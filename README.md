# sefa
## Simulation
* [Implement a faulty AES For SIFA and SEFA](https://github.com/sefasifa/sefa/blob/main/README.md#implement-a-faulty-aes-for-sifa-and-sefa) 
  * [Regular Fault](https://github.com/sefasifa/sefa/blob/main/README.md#regular-fault) 
  * [Missrate](https://github.com/sefasifa/sefa#Missrate)
  * [Dummy Round](https://github.com/sefasifa/sefa#dummy-round)
  * [Error-Correction mode](https://github.com/sefasifa/sefa#error-correction-mode) 
* [SEI and LLR Computation](https://github.com/sefasifa/sefa#sei-and-llr-computation)
* [Key-Recovery for a byte](https://github.com/sefasifa/sefa#key-recovery) 


**Random Fault Value**

 First of all we define fault value. Here we used the random number generator which is provided by matlab for random fault,. Since the Fault which is propsed in the attack is 
(fb) bits *random And*. 

```matlab
faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
```
### Implement a faulty AES For SIFA and SEFA 

[AES code](https://nevonprojects.com/aes-source-code-inmatlab/)  is written by J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

```matlab
function [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex)
%AES_INIT  Initialisation of AES-components.


% Create the S-box and the inverse S-box
[s_box, inv_s_box] = s_box_gen (1);

% Create the round constant array
rcon = rcon_gen (1);

% % Convert the cipher key from hexadecimal (string) to decimal representation
   key = hex2dec(key_hex);


% Create the expanded key (schedule)
w = key_expansion (key, s_box, rcon, 1);

% Create the polynomial transformation matrix and the inverse polynomial matrix
% to be used in MIX_COLUMNS
[poly_mat, inv_poly_mat] = poly_mat_gen (1);
end

function state_out = add_round_key (state_in, round_key)
%ADD_ROUND_KEY  Add (XOR) the round key to the state.
% Add state (matrix) and round key (matrix) via bitwise XOR
state_out = bitxor (state_in, round_key);


% Apply the just created index matrix to the input matrix.
% Elements of the index matrix are linear (column-wise) indices.
matrix_out = matrix_in (ind_mat);
function state_out = shift_rows (state_in)
%SHIFT_ROWS  Cyclically shift the rows of the state matrix.

function state_out = mix_columns (state_in, poly_mat)
%MIX_COLUMNS  Transform each column of the state matrix.
% Define the irreducible polynomial 
% to be used in the modulo operation in poly_mult
mod_pol = bin2dec ('100011011');

% Loop over all columns of the state matrix
for i_col_state = 1 : 4
        
    % Loop over all rows of the state matrix
    for i_row_state = 1 : 4

        % Initialize the scalar product accumulator
        temp_state = 0;
        
        % For the (innner) matrix vector product we want to do
        % a scalar product 
        % of the current row vector of poly_mat
        % and the current column vector of the state matrix.
        % Therefore we need a counter over 
        % all elements of the current row vector of poly_mat and
        % all elements of the current column vector of the state matrix
        for i_inner = 1 : 4
        
            % Multiply (GF(2^8) polynomial multiplication)
            % the current element of the current row vector of poly_mat with
            % the current element of the current column vector of the state matrix
            temp_prod = poly_mult (...
                        poly_mat(i_row_state, i_inner), ...
                        state_in(i_inner, i_col_state), ...
                        mod_pol);
            
            % Add (XOR) the recently calculated product
            % to the scalar product accumulator
            temp_state = bitxor (temp_state, temp_prod);
                        
        end
        
        % Declare (save and return) the final scalar product accumulator
        % as the current state matrix element
        state_out(i_row_state, i_col_state) = temp_state;
        
    end
    
end

function [s_box, inv_s_box] = s_box_gen (vargin)

if nargin > 0
    
    % Switch the verbose mode flag on
    verbose_mode = 1;
    
% If there is no optional "verbose mode" argument
else
    
    % Switch the verbose mode flag off
    verbose_mode = 0;
    
end


mod_pol = bin2dec ('100011011');

% The polynomial multiplicative inverse of zero is defined here as zero.
% Matlab vectors start with an index of "1"
inverse(1) = 0;
load('s_box.mat');

inv_s_box = s_box_inversion (s_box);

% Display intermediate result if requested
if verbose_mode
    

    s_box_mat = reshape (s_box, 16, 16)';
%     disp_hex ('    s_box : ', s_box_mat)
    inv_s_box_mat = reshape (inv_s_box, 16, 16)';
%     disp_hex ('inv_s_box : ', inv_s_box_mat)
    
end
```
**Faulty & Non-Faulty Cipher**

Faults are injected at the beginning of round 10 in these three situations: Regular fault, missrate faults and Error-Correcting mode. In Dummy-Rounds, 
faults can be injected in different rounds, or even in dummies. In cipher Functions, *faultee* shows faulty or non-faulty computations. When *faultee==0*, the cipher is correct
and fault is not injected.[Cipher](https://nevonprojects.com/aes-source-code-inmatlab/) 

```matlab
function [ciphertext,partialstate,round_key] = cipher (plaintext, w, s_box, poly_mat, vargin,faultee,faultvalue)

partialstate=zeros(4,4,9);
%   Version 1.0     30.05.2001

% If there is an optional "verbose mode" argument
if nargin > 4
    
    % Switch the verbose mode flag on
    verbose_mode = 1;
    
% If there is no optional "verbose mode" argument
else
    
    % Switch the verbose mode flag off
    verbose_mode = 0;
    
end

% If the input vector is a cell array or does not have 16 elements
if iscell (plaintext) | prod (size (plaintext)) ~= 16

    % Inform user and abort
    error ('Plaintext has to be a vector (not a cell array) with 16 elements.')
    
end

% If any element of the input vector cannot be represented by 8 bits
if any (plaintext < 0 | plaintext > 255)
    
    % Inform user and abort
    error ('Elements of plaintext vector have to be bytes (0 <= plaintext(i) <= 255).')
    
end

% If the expanded key array is a cell arrray or does not have the correct size
if iscell (w) | any (size (w) ~= [44, 4])

    % Inform user and abort
    error ('w has to be an array (not a cell array) with [44 x 4] elements.')
    
end

% If any element of the expanded key array can not be represented by 8 bits
if any (w < 0 | w > 255)
    
    % Inform user and abort
    error ('Elements of key array w have to be bytes (0 <= w(i,j) <= 255).')
    
end


state = reshape (plaintext, 4, 4);


round_key = (w(1:4, :))';


state = add_round_key (state, round_key);

% Loop over 9 rounds
for i_round = 1 : 9

     state = sub_bytes (state, s_box);

    state = shift_rows (state);
    


    state = mix_columns (state, poly_mat);

    round_key = (w((1:4) + 4*i_round, :))';

%% Faults are injected here for  Dummy-Round.

    if (faultee-10)==i_round
         state(1,1)=bitand(state(1,1),faultvalue);
    end
    state = add_round_key (state, round_key);

end
%% Faults are injected here for  Regular fault, missrate faults and Error-Correcting mode.
    if i_round==9
        if (faultee>0)
            state(1,1)=bitand(state(1,1),faultvalue);
        end;
    end
state = sub_bytes (state, s_box);


state = shift_rows (state);
    

round_key = (w(41:44, :))';

state = add_round_key (state, round_key);
    
ciphertext = reshape (state, 1, 16);

```
#### Regular-Fault 
In all cases, key and palintext is considered as random. The total number of random key which is considered is key_t and  sample_t show number of random plaintext for each key_t. For each sample, faulty *(cipherc)* and non-faulty *(cipherf)* encryption is computed. In any faulty encryption *fb* bits is injected.   

```matlab
function [key_col,cipherc,cipherf]=regularfault(sample_t,key_t,fb)
key_col=cell(1,1000);
cipherc=zeros(16,15000,100);
cipherf=zeros(16,15000,100);
    for key_n=1:key_t
        key(1,1)=round(rand*255)  ;
        key(2,1)=round(rand*255)  ;
        key(3,1)=round(rand*255)  ;
        key(4,1)=round(rand*255)  ;
        key(5,1)=round(rand*255)  ;
        key(6,1)=round(rand*255)  ;
        key(7,1)=round(rand*255)  ;
        key(8,1)=round(rand*255)  ;
        key(9,1)=round(rand*255)  ;
        key(10,1)=round(rand*255) ;
        key(11,1)=round(rand*255) ;
        key(12,1)=round(rand*255) ;
        key(13,1)=round(rand*255) ;
        key(14,1)=round(rand*255) ;
        key(15,1)=round(rand*255) ;
        key(16,1)=round(rand*255) ;
        key_col{key_n} = {dec2hex(key)};
        SS=reshape(key,1,16);
        key_hex=dec2hex(SS);
        [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex);
            for i=1:sample_t
                faultee=1;
                plaintext(1,1)=round(rand*255);
                plaintext(2,1)=round(rand*255);
                plaintext(3,1)=round(rand*255);
                plaintext(4,1)=round(rand*255);
                plaintext(5,1)=round(rand*255);
                plaintext(6,1)=round(rand*255);
                plaintext(7,1)=round(rand*255);
                plaintext(8,1)=round(rand*255);
                plaintext(9,1)=round(rand*255);
                plaintext(10,1)=round(rand*255);
                plaintext(11,1)=round(rand*255);
                plaintext(12,1)=round(rand*255);
                plaintext(13,1)=round(rand*255);
                plaintext(14,1)=round(rand*255);
                plaintext(15,1)=round(rand*255);
                plaintext(16,1)=round(rand*255);
                faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
                [ciphertext,] = cipher (plaintext, w, s_box, poly_mat,0,0,faultvalue);
                cipherc(:,i,key_n)=ciphertext;
                [ciphertextf,] = cipher (plaintext, w, s_box, poly_mat,0,faultee,faultvalue);
                cipherf(:,i,key_n)=ciphertextf;
            end
    end
end
```


#### Missrate 
In all cases, key and palintext is considered as random. We define missrate such that, a fault can not be injected. So in this case, faultee should be 0. We define a random function to create random 0 or 1 value for desired missrate.

```matlab
Pe=100*missrate;
random_f= round(Pe*rand(key_t,sample_t));
random_f(:)=(random_f(:)<100-Pe);
exactmissrate=sum(sum(random_f(1:key_t,:)))/sample_t*key_t;
```
#### Dummy Round 
 In Dummy rounds, we create a random vector by K*10 elements. Then 10 random round is chosen and sorted. If the considered SEL_R which we injected faults is equal by the last sorted chosen round, then fault is useful; otherwise the fault is injected to other rounds or in dummy rounds. If fault is injected in dummy, *faultee* shoulde be equal to zero, otherwise the round of injecteion is added to *faultee* to specifty the number of round in **cipher** function. 

```matlab
  selected_rand=sort(round(randperm(k*10,10)));
  if all(selected_rand(:)~=sel_R)
      faultee=0;
  elseif selected_rand(10)==sel_R
      faultee=1;
  else
      fault_R=find(selected_rand==sel_R);
      faultee=fault_R+10;
  end
```
#### Error-Correction Mode
We consider that in this mode d' bit faults can be corrected. In this case, when the HW of injection is less than d', plaintexte encrypted in non-faulty mode. 

```matlab
    hw=8-sum(data_fault(:,:)==1);
    if hw<=d'
        faultee=0;
    end
```



### SEI and LLR Computation
 In this part we compute the inverse of CIPHERC to the input of Sbox in the begining of the round 10 by guessing 256 key-guesses. Then the computed inverse would be masked by AND. 
```matlab
for key_g=0:255
                if isequal(cipherc(:,i,key_n),cipherf(:,i,key_n))
                    data_d_coll_i(key_n,key_g+1,ineff)=bitand(sub_bytes(bitxor((cipherc(1,i,key_n)),key_g),inv_s_box),mask);
                    delta_d_i(key_n,data_d_coll_i(key_n,key_g+1,ineff)+1,key_g+1)=delta_d_i(key_n,data_d_coll_i(key_n,key_g+1,ineff)+1,key_g+1)+1;
                else
                    data_d_coll_e(key_n,key_g+1,eff)=bitand(sub_bytes(bitxor((cipherc(1,i,key_n)),key_g),inv_s_box),mask);
                    delta_d_e(key_n,data_d_coll_e(key_n,key_g+1,eff)+1,key_g+1)=delta_d_e(key_n,data_d_coll_e(key_n,key_g+1,eff)+1,key_g+1)+1;
                end 
            end 
            for m=1:(mask+1)
                if ni>0
                    SEI_i(key_n,i,:)=(((delta_d_i(key_n,m,:)/ni)-(1/(mask+1))).^2)+SEI_i(key_n,i,:);
                end
                if ne>0
                    SEI_e(key_n,i,:)=(((delta_d_e(key_n,m,:)/ne)-(1/(mask+1))).^2)+SEI_e(key_n,i,:);
                    if ni>0
                        for ms=1:(mask+1)
                            SEI_joint(key_n,i,:)=((((delta_d_i(key_n,ms,:)).*delta_d_e(key_n,m,:)/(ne*ni))-(1/((mask+1)^2))).^2)+SEI_joint(key_n,i,:);
                        end
                    end
                end
                for key_g=1:256
                     if (delta_d_e(key_n,m,key_g))>0 && ne>0
                        LLR_E(key_n,i,key_g)=ne*q_e(m)*log2((delta_d_e(key_n,m,key_g)/ne)/(1/(mask+1)))+ LLR_E(key_n,i,key_g);
                     end
                     if (delta_d_i(key_n,m,key_g))>0 && ni>0
                        LLR_I(key_n,i,key_g)=ni*q_i(m)*log2((delta_d_i(key_n,m,key_g)/ni)/(1/(mask+1)))+LLR_I(key_n,i,key_g);
                        if (delta_d_e(key_n,m,key_g))>0 && ne>0
                            for ms=1:(mask+1)
                                            LLR_joint(key_n,i,key_g)=ni*q_i(ms)*ne*q_e(m)*log2((((delta_d_i(key_n,ms,key_g)).*delta_d_e(key_n,m,key_g)/(ne*ni))/(1/((mask+1)^2))))+LLR_joint(key_n,i,key_g);
                            end
                        end
                     end
                end
            end

```
### Key Recovery
```matlab
   rank_eff_sei(key_n,i)=sum(SEI_e(key_n,i,:)>=SEI_e(key_n,i,key_b_0+1));
   rank_ineff_sei(key_n,i)=sum(SEI_i(key_n,i,:)>=SEI_i(key_n,i,key_b_0+1));
   rank_joint_sei(key_n,i)=sum(SEI_joint(key_n,i,:)>=SEI_joint(key_n,i,key_b_0+1));
   rank_eff_llr(key_n,i)=sum(LLR_E(key_n,i,:)>=LLR_E(key_n,i,key_b_0+1));
   rank_ineff_llr(key_n,i)=sum(LLR_I(key_n,i,:)>=LLR_I(key_n,i,key_b_0+1)); 
   rank_joint_llr(key_n,i)=sum(LLR_joint(key_n,i,:)>=LLR_joint(key_n,i,key_b_0+1));
```




