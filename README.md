# Processing WAV files with Fast Fourier Transforms in C++
**Jennifer Freeman**

## Description
This program implements the Cooley-Tukey radix-2 fast Fourier Transform (FFT) on a WAV audio file. The user can apply audio processing effects to output a modified WAV audio file.
### Cooley-Tukey radix-2 FFT algorithm

An FFT is as algorithm for efficiently computing the Discrete Fourier Transform (DFT) on a data set. For a vector of data $x$ of size $N$ the DFT is defined as:
$$X[k]=\sum_{n=0}^{N-1} x[n]e^{\frac{-2\pi i nk}{N}}$$
The Cooley-Tukey radix-2 FFT algorithm computes the FFT for a complex data set of size $N=2^m$. Given this N we can partition the DFT into two parts.

$$X[k]=\underbrace{\sum_{m=0}^{N/2-1} x[2m]e^{\frac{-2\pi i mk}{N/2}}}_{\text{$x_{even}$}} + \underbrace{e^{\frac{-2\pi i k}{N}}}_{\text{$W_N^k$}}\underbrace{\sum_{m=0}^{N/2-1} x[2m+1]e^{\frac{-2\pi i mk}{N/2}}}_{\text{$x_{odd}$}}$$

The first summation $x_{even}$ contains the even indices of $x$, $x_{odd}$ contains the odd indices, and $W_N^k$ is called the *twiddle* factor. Using the fact that $e^{z+2\pi i} = e^z$, then the following relationship holds:

$$X[k] = x_{even}^k + W_N^k x_{odd}^k, \quad X[k+N/2] = x_{even}^k - W_N^k x_{odd}^k$$

This saves additional computations since each $x_{even}^k$ $x_{odd}^k$ can be used for two different DFT calculations. Since $N$ is a power of 2, $N$ can be successively halved and the algorithm recursively computes the DFT [^1]. This recursive quality leads to a computationally efficient algorithm.

### Digital filters

A digital filter is an operation in signal processing that alters specific frequency components within a signal [^2]. In an audio signal, we can apply digital filters to reduce or enhance specific frequencies.

A low-pass digital filter attenuates frequencies above some threshold frequency and allows all lower frequencies to "pass" through the filter. The simplest low-pass frequency is called a *sharp cut-off* filter where all frequencies above the threshold are removed and all other frequencies remain unchanged [^2].

To apply the low-pass sharp cut-off filter to an audio file, an FFT is performed to produce the frequency spectrum $X[k]$ of size $N$. A frequency $f$ occurs at the $k$-th value of $X$ and satisfies $f = f_{samp}k/N$ where $f_{samp}$ is the sampling rate of the the audio signal. For a given frequency, $k$ can be solved for and the resulting $X[k]$ values can be set to zero. The inverse FFT is then performed on the modified $X[k]$ and the resulting filtered audio signal is obtained. 


## Program Input

The program takes three arguments on input in the following order:
1. A WAV audio file name with a maximum file size of $262144*2+44=524332$ bytes. Example: "hello_world.wav"
2. The digital filter type to apply to the audio file. Valid inputs are: "low"
3. The frequency in Hertz to use with the filter between the audible range of 20 to 20000. Example: 2400

## Program Output

The program outputs a WAV audio file with the transformed audio information from the input file. The output file is named with the system data and time. 

## References

[^1]: Bekele, A. J. Advanced Algorithms. (2016). Cooley-tukey fft algorithms. Advanced algorithms.

[^2]: O'Haver, T. (2019, December). Fourier filter. Retrieved 16:00, December 12, 2021, from https://terpconnect.umd.edu/~toh/spectrum/FourierFilter.html