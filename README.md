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
%
%   STATE_OUT = MIX_COLUMNS (STATE_IN, POLY_MAT) 
%   operates on the state matrix STATE_IN column-by-column
%   using POLY_MAT as the transformation matrix.
%
%   MIX_COLUMNS can also directly compute 
%   the inverse column transformation INV_MIX_COLUMNS
%   by utilizising the inverse transformation matrix INV_POLY_MAT.

%   Copyright 2001-2005, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.0     30.05.2001

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












state_out = cycle (state_in, 'left');



```



