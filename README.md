# Schwarzschild QNMs via a Leaver-type Series Method
Here I will use a Frobenius-series (Leaver-type) solver for Schwarzschild black hole (BH) quasinormal modes (QNMs), which will be Step 1 of a project that aims to compare this numerical method to the semiclassical CFT/Nekrasov-Shatashvili quantization condition (this will form part of Step 2). 
## Python Scripts
- * `schwarzschild.py` : The solver. In this script, the Frobenius recursion relation is derived for the Regge-Wheeler equation with spin s=0,1,2 and angular momentum number l. A truncated, row-normalized matrix M is constructed and the QNM frequencies are found from the roots of its determinant. Specifically, a Hill determinant method was used, similar to Leaver's continued-fraction method.
