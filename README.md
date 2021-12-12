# Processing WAV files with Fast Fourier transform in C++
**Jennifer Freeman**

## Description
This program implements the Cooley-Tukey radix-2 fast Fourier transform (FFT) on a WAV audio file.

## Cooley-Tukey radix-2 FFT algorithm

An FFT is as algorithm for efficiently computing the Discrete Fourier Transform (DFT) on a data set. For a vector of data $x$ of size $N$ the DFT is defined as:
$$X[k]=\sum_{n=0}^{N-1} x[n]e^{\frac{-2\pi i nk}{N}}$$
The Cooley-Tukey radix-2 FFT algorithm computes the FFT for a complex data set of size $N=2^m$. Given this N we can partition the DFT into two,

$$X[k]=\underbrace{\sum_{m=0}^{N/2-1} x[2m]e^{\frac{-2\pi i mk}{N/2}}}_{\text{$x_{even}$}} + \underbrace{e^{\frac{-2\pi i k}{N}}}_{\text{$W_N^k$}}\underbrace{\sum_{m=0}^{N/2-1} x[2m+1]e^{\frac{-2\pi i mk}{N/2}}}_{\text{$x_{odd}$}}$$

The first summation $x_{even}$ contains the even indices of $x$, $x_odd$ contains the odd indices, and $W_N^k$ is called the *twiddle* factor.



## Program Input

The program takes on input a WAV audio file. 
--- padding

## Program Output


[^1]: Bekele, A. J. A. A. (2016). Cooley-tukey fft algorithms. Advanced algorithms.