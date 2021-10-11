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
end


% Apply the just created index matrix to the input matrix.
% Elements of the index matrix are linear (column-wise) indices.
matrix_out = matrix_in (ind_mat);
end
function state_out = shift_rows (state_in)
%SHIFT_ROWS  Cyclically shift the rows of the state matrix.
%

state_out = cycle (state_in, 'left');



```



