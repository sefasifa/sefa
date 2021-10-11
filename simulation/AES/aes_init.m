function [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex)
%AES_INIT  Initialisation of AES-components.
%
%   [S_BOX, INV_S_BOX, W, POLY_MAT, INV_POLY_MAT] = AES_INIT
%   initializes AES-components 
%   to be used by subsequent functions.
%
%   In the initialization step the S-boxes (S_BOX and INV_S_BOX) 
%   and the polynomial matrices (POLY_MAT and INV_POLY_MAT)
%   are created and
%   an example cipher key is expanded into 
%   the round key schedule (W).

%   Copyright 2001-2005, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.0     30.05.2001

% Clear the command window
% clc

% Create the S-box and the inverse S-box
[s_box, inv_s_box] = s_box_gen (1);

% Create the round constant array
rcon = rcon_gen (1);

% Define an arbitrary 16-byte cipher key in hexadecimal (string) representation
% The following two specific keys are used as examples 
% in the AES-Specification (draft)
% key_hex = {'11' '01' '02' '23' '04' '05' '06' '07' ...
%           '08' '09' '0a' '0b' '0c' '0d' '0e' '0f'};
% key_hex = {'2b' '7e' '15' '16' '28' 'ae' 'd2' 'a6' ...
%            'ab' 'f7' '15' '88' '09' 'cf' '4f' '3c'};
  
% key_hex = {'00' '00' '00' '00' '00' '00' '00' '00' ...
%            '00' '00' '00' '00' '00' '00' '00' '0d'};
% 
%  key_hex = {'56' '3c' '37' '27' 'b1' '4a' '75' 'ff' ...
%             '56' 'd4' '67' 'ae' '5d' 'da' '6f' 'fe'};    
% % 
%  key_hex = {'58' '34' '3c' '27' 'b1' '4a' '75' 'ff' ...
%             '57' 'd6' '65' 'ac' 'ad' 'da' '6f' 'fe'};    
% % Convert the cipher key from hexadecimal (string) to decimal representation
   key = hex2dec(key_hex);
%     key(1,1)=round(rand*255)  ;
%     key(2,1)=round(rand*255)  ;
%     key(3,1)=round(rand*255)  ;
%     key(4,1)=round(rand*255)  ;
%     key(5,1)=round(rand*255)  ;
%     key(6,1)=round(rand*255)  ;
%     key(7,1)=round(rand*255)  ;
%     key(8,1)=round(rand*255)  ;
%     key(9,1)=round(rand*255)  ;
%     key(10,1)=round(rand*255) ;
%     key(11,1)=round(rand*255) ;
%     key(12,1)=round(rand*255) ;
%     key(13,1)=round(rand*255) ;
%     key(14,1)=round(rand*255) ;
%     key(15,1)=round(rand*255) ;
%     key(16,1)=round(rand*255) ;
%     key(1,1)=round(rand*255)  
%     key(2,1)=round(rand*255)  
%     key(3,1)=round(rand*255)  
%     key(4,1)=round(rand*255)  
%     key(5,1)=round(rand*255)  
%     key(6,1)=round(rand*255)  
%     key(7,1)=round(rand*255)  
%     key(8,1)=round(rand*255)  
%     key(9,1)=round(rand*255)  
%     key(10,1)=round(rand*255) 
%     key(11,1)=round(rand*255) 
%     key(12,1)=round(rand*255) 
%     key(13,1)=round(rand*255) 
%     key(14,1)=round(rand*255) 
%     key(15,1)=round(rand*255) 
%     key(16,1)=round(rand*255) 
%     key(1,1) = 247 ;
%     key(2,1) =  70 ;
%     key(3,1) =  63 ;
%     key(4,1) = 171 ;
%     key(5,1) =  68 ;
%     key(6,1) =  69 ;
%     key(7,1) = 128 ;
%     key(8,1) = 166 ;
%     key(9,1) =  83 ;
%     key(10,1)=  46 ;
%     key(11,1)= 194 ;
%     key(12,1)= 242 ;
%     key(13,1)= 234 ;
%     key(14,1)= 215 ;
%     key(15,1)= 219 ;
%     key(16,1)= 243 ;
%     key(1,1) =  237  ;
%     key(2,1) =  180  ;
%     key(3,1) =  193  ;
%     key(4,1) =  247  ;
%     key(5,1) =  122  ;
%     key(6,1) =   31  ;
%     key(7,1) =  163  ;
%     key(8,1) =  217  ;
%     key(9,1) =   58  ;
%     key(10,1)=  173  ;
%     key(11,1)=  254  ;
%     key(12,1)=   29  ;
%     key(13,1)=   65  ;
%     key(14,1)=  243  ;
%     key(15,1)=  197  ;
%     key(16,1)=    4  ;

% Create the expanded key (schedule)
w = key_expansion (key, s_box, rcon, 1);

% Create the polynomial transformation matrix and the inverse polynomial matrix
% to be used in MIX_COLUMNS
[poly_mat, inv_poly_mat] = poly_mat_gen (1);