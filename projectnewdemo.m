clear all;
close all;
clc

%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%
L = 56; %datapacket bitsize
P = 300; % datapacket amount must be even
startdata = randi([0 1],L,P);

N0 = 0.4;
Es = 1;
variance = N0/2;
%%%%%%%Hamming encoding%%%%%%%
hamming = hammingencoding(startdata);

%%%%%%%interleaving%%%%%%%
interleaved = interleaving(hamming);
lengthofdata = length(interleaved);
X = sprintf('Data packet final length: %d', lengthofdata);
disp(X)
payload = P*L;
Y = sprintf('Data packet payload length: %d', payload);
disp(Y)
overhead = lengthofdata-payload;
overheadprocent = round((overhead)/lengthofdata*100);
Z = sprintf('Overhead payload size: %i bits, which is around %i%%', overhead, overheadprocent);
disp(Z)

%%%%%%%%%%%vectorencode%%%%%%%%%%%
vector = vectorencode(interleaved);
normalized = normalize(vector,Es);
%%%%%%%%%%%Channel%%%%%%%%%%%%
noised = addnoise(normalized,variance);
%plotthis = plotting(noised);
%%%%%%%%%%after channel%%%%%%%%%%%%
%%%%%%%%%%vector decode%%%%%%%%%%%%
received = vectordecode(noised);
%%%%%%%%%%deinterleave%%%%%%%%%%%%
deinterleaved = deinterleave(received,P);
%%%%%%%%%%dehamming%%%%%%%%%%%%
dehamming = hammingdecode(deinterleaved);

fejludenhammingdecode = checkfejl(hamming,deinterleaved);
fejlprocentudenhamming = fejludenhammingdecode/lengthofdata*100
fejlmedhamming = checkfejl(startdata,dehamming);
fejlprocentmedhamming = fejlmedhamming/(L*P)*100


%%%%%%%%%%%%%%%%%%%%%%FUCKTIONER%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%skal stå til sidst i dokumentet%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% create matrix %%%%%%%%%%%%%%%%%%
function matx = getmatrix(x,L)
   z = zeros(x,L+x);
   d = (1:L+x);
   b = de2bi(d);
   matx = transpose(b);

end

%%%%%%%%%%%Hamming function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hamming = hammingencoding(data)
      sizes = size(data);
      L = sizes(1);
      P = sizes(2);
      minimumr = log2(L+log2(L+1)+1);
      PBits = ceil(minimumr);
      packetlength = L+PBits;
      potens = 0;
      for j=1:1:P
         for i=1:1:packetlength
             if i == 2^potens
                     hamming(i,j) = 2;
                     potens = potens + 1;
               else
                    hamming(i,j) = data(i-potens,j);
             end
         end
         potens = 0;
         %%%%%%%parity check%%%%%%%%
         matx = getmatrix(PBits,L);
         pcheck = matx * hamming(:,j);
         pos = getparitypos(PBits);
         hamming(:,j) = insertparity(pos,PBits,pcheck,hamming(:,j));
      end


end
%%%%%%%%%%%%%%%%%%%Insert parity values%%%%%%%%%%%%%
function hamming = insertparity(pos,PBits,pcheck,hamming)
   for i = 1 :1 : PBits
   value = mod(pcheck(i,1),2);
   hamming(pos(i),1) = value;
   end
end

%%%%%%%%%%% Get parity positions%%%%%%%%%%%%
function pos = getparitypos(PBits)
  for i=1:1:PBits
    pos(i,1) = 2^(i-1);
  end
end

%%%%%%%%%%%%%%%%%%%%Interleaving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interleaved = interleaving(data)
   sizes = size(data);
   L = sizes(1);
   P = sizes(2);
   interleaved = [];
   for i=1:1:L
       for j=1:1:P
        interleaved(end+1,1) = data(i,j);
       end
   end
end
%%%%%%%%%%%%%%%%%%%Vector Encoding%%%%%%%%%%%%%%%%%%%%%%
function vector = vectorencode(data)
   L = length(data);
   for i=1:2:L
       val = data(i,1);
       val2 = data(i+1,1);

       switch val
           case 1
               switch val2
                   case 1
                       %if data == 11 vector == 1,1
                       vector(i,1) = 1;
                       vector(i+1,1) = 1;
                   case 0
                       %if data == 10 vector == 1,-1
                       vector(i,1) = 1;
                       vector(i+1,1) = -1;
               end
           case 0
               switch val2
                   case 1
                       %if startdata == 01 vector == -1,1
                       vector(i,1) = -1;
                       vector(i+1,1) = 1;
                   case 0
                       %if startdata == 00 vector == -1,-1
                       vector(i,1) = -1;
                       vector(i+1,1) = -1;
               end
       end

   end
end
%%%%%%%%%%%%%%%%%%%%%%Normalize%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normalized = normalize(vector,Es)
   normalized = (Es)^0.5 *vector;
end

%%%%%%%%%%%%%%%%%%%%Add Noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noised = addnoise(normalized,var)
      lengthofdata = length(normalized);
      tal = randn(1,lengthofdata);
      noise = var^0.5*tal;
      noised = normalized + transpose(noise);
      %histogram(noise(1,1:lengthofdata),20);
end
%%%%%%%%%%%%%%%%%%%%plot Noise received %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotthis = plotting(data)
   L = length(data);
   val = [];
   val2 = [];
   for i=1:2:L
       val(end+1,1) = data(i,1);
       val2(end+1,1) = data(i+1,1);

   end
   plot(val,val2,".")
   plotthis = 0;
end
%%%%%%%%%%%%%%%%%%%% Vector decoding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function received = vectordecode(data);
   L = length(data);
   for i=1:2:L
       val = data(i,1);
       val2 = data(i+1,1);
       if val >= 0
            if val2 >= 0
               received(i,1) = 1;
               received(i+1,1) = 1;
            elseif val2 < 0
               received(i,1) = 1;
               received(i+1,1) = 0;
            end

       elseif val < 0
            if val2 >= 0
              received(i,1) = 0;
              received(i+1,1) = 1;
            elseif val2 < 0
              received(i,1) = 0;
              received(i+1,1) = 0;
            end
        end


   end
end
%%%%%%%%%%%%%%%%%%%%Deinterleaving%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function deinter = deinterleave(data,num)
    L = length(data)/num;
    P = length(data)/L;
    i = 1;
        for k=1:1:L
            for j=1:1:P
                    deinter(k,j) = data(i);
                    i = i+1;
            end
        end
 end
 %%%%%%%%%%%%%%%%%%%%Dehamming%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dehamming = hammingdecode(data)
      corrected = errorcorrection(data);
      sizes = size(data);
      L = sizes(1);
      P = sizes(2);
      minimumr = log2(L);
      PBits = ceil(minimumr);
      pos = getparitypos(PBits);
      dehamming = removePBits(corrected,pos);
end

%%%%%%%%%%%%%%%%%%%%Correction of errors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrected = errorcorrection(data)
   sizes = size(data);
   L = sizes(1);
   P = sizes(2);
   minimumr = log2(L);
   PBits = ceil(minimumr);
   matx = getmatrix(PBits,L-PBits);
   for i=1:1:P
        tempdata = data(:,i);
        pcheck = transpose(mod(matx * tempdata,2));
        errorpos = bi2de(pcheck);
       if errorpos ~= 0 && errorpos <= L
           if tempdata(errorpos,1) == 1
               tempdata(errorpos,1) = 0;
            elseif tempdata(errorpos,1) == 0
              tempdata(errorpos,1) = 1;
           end
       end
       corrected(:,i) = tempdata;
   end
end


function dehamming = removePBits(data, pos)
   sizes = size(data);
   L = sizes(1);
   P = sizes(2);
   amount = length(pos);
   for j=1:1:P
       tempdata = data(:,j);
     for i=1:1:amount
        position = pos(i,1);
        tempdata(position)=3;
     end
     tempdata = tempdata(tempdata~=3);
     dehamming(:,j)=tempdata;
   end
end

%%%%%%%%%%%Fejl tjekker%%%%%%%%%%
function fejl = checkfejl(data,data2)
   mat = data-data2;
   wups = size(mat);
   n=0;
   for i=1:1:wups(1)
      for j=1:1:wups(2)
        if mat(i,j) ~= 0
          n = n+1;
        end
      end
       end
   fejl = n;
end
