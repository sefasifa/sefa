# sefa
## Simulation
* Implement a faulty AES For SIFA and SEFA 
  * Regular Fault 
  * Fault by considering Missrate
  * Fault by considering Dummy Round
  * Fault in Error-Correction mode 
* SEI and LLR calculation
* Key-Recovery for a byte 


**Random Fault Value**

 First of all we define fault value. Here we used the random number generator which is provided by matlab for random fault,. Since the Fault which is propsed in the attack is 
(fb) bits *random And*. 

```matlab
faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
```
### Implement a faulty AES For SIFA and SEFA 

AES code  is written by J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

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
**This part is used to creates a faulty/non-faulty cipher**
Faults are injected at the beginning of round 10 in these three situations: Regular fault, missrate faults and Error-Correcting mode. In Dummy-Rounds situation, 
faults can be injected in different rounds, or even in dummies. In cipher Functions, *faultee* shows faulty or non-faulty computations. 

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













