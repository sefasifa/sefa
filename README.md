# sefa
## Simulation
* Implement a faulty AES For SIFA and SEFA 
  * Regular Fault 
  * Fault by considering Missrate
  * Fault by considering Dummy Round
  * Fault in Error-Correction mode 
* SEI and LLR calculation
* Key-Recovery for a byte 


**Implement a faulty AES For SIFA and SEFA **

For random fault, here we used the random number generator which is provided by matlab. Since the Fault which is propsed in the attack is 
(fb) bits *random And*, the definition based on the random fault is 

```matlab
faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
'''



