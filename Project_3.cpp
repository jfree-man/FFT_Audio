/**
 * @file main.cpp
 * @author Jennifer Freeman (freemjc@mcmaster.ca)
 * @brief 
 * @version 0.1
 * @date 2021-12-12
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <complex>
#include <chrono>
#include <ctime>
#include <iomanip>
using namespace std;

/**
 * @brief WAV audio file class.
 * 
 */
class wav_file
{
public:
    // ===========
    // Constructor
    // ===========

    /**
      * @brief Construct a new wav file object.
      * 
      * @param input_string string of wav file name, (ex. "filename.wav").
      */
    wav_file(const string &input_string)
    {
        read_wav(input_string);
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Read in WAV file.
     * 
     * @param input_string string of wav file name, (ex. "filename.wav").
     */
    void read_wav(const string &input_string)
    {

        string s;                        // String to store input.
        string wav_identifier;           // String to store expected WAV identifier.
        string current;                  // Current input value.
        int16_t temp;                    // Temporary variable to store input values.
        uint32_t sampling_rate;          // Sampling rate in WAV file.
        uint64_t data_section = 0;       // Size of data section in file.
        uint64_t power = 1;              // Nearest power of 2 to data size.
        uint64_t data_size = 1;          // Size of data in 2-bytes.
        uint64_t max_data_size = 262144; // Maximum allowable data size (2^18).

        // Open input file.
        ifstream input(input_string, ios::binary);

        if (!input.is_open())
        {
            throw ios_base::failure("Error opening input file! \n");
        };

        // Seek to WAV identifier position.
        input.seekg(8);

        // Get WAV identifier.
        for (uint32_t i = 0; i < 4; i++)
        {
            current = (char)input.get();
            wav_identifier.insert(i, current);
        }

        if (wav_identifier != "WAVE")
        {
            throw invalid_argument("Not a .WAV file.");
        }

        // Read in WAV header (44 bytes).
        input.seekg(0);
        getline(input, s);
        header = s.substr(0, 44);

        // Get data section size.
        data_section = *(int32_t *)&header[40];

        // Ensure data section is a multiple of 2.
        if (data_section % 2 != 0)
        {
            data_section = data_section + 1;
        }

        // Calculate size of data in 16 bits.
        data_size = data_section / 2;

        if (data_size > max_data_size)
        {
            throw invalid_argument("File size is larger than maximum allowable.");
        }

        // Compute next closest power of two to data size.
        while (power <= data_size)
        {
            power *= 2;
        }

        // Seek to data chunk.
        input.seekg(44);

        // Read in data.
        for (uint64_t i = 0; i < data_size; i++)
        {
            input.read((char *)&temp, sizeof(int16_t));
            data.push_back((int16_t)temp);
        }

        // Resize data to include padding.
        data.resize(power);

        input.close();

        // Get sampling rate.
        sampling_rate = *(uint32_t *)&header[24];
        f_samp = sampling_rate;
    }

    /**
     * @brief Size of data in WAV file.
     * 
     * @return uint64_t
     */
    uint64_t datasize()
    {
        return data.size();
    }

    /**
     * @brief Get the data object and convert to complex<double>.
     * 
     * @return vector<complex<double>> 
     */
    vector<complex<double>> get_data()
    {
        vector<complex<double>> new_vec(data.begin(), data.end());
        return new_vec;
    }

    /**
     * @brief Get the header object.
     * 
     * @return string 
     */
    string get_header()
    {
        return header;
    }

    /**
     * @brief Get the sampling rate object.
     * 
     * @return uint32_t 
     */
    uint32_t get_samp()
    {
        return f_samp;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Header in WAV file.
     * 
     */
    string header;

    /**
     * @brief Data in WAV file.
     * 
     */
    vector<int16_t> data;

    /**
     * @brief sampling rate
     * 
     */
    uint32_t f_samp;
};

/**
 * @brief Fast Fourier transform class.
 * 
 */
class fft
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new fft object.
     * 
     * @param dat data to perform FFT on.
     */
    fft(vector<complex<double>> &dat)
    {
        x = cooley_tukey(dat);
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Cooley-Tukey radix-2 FFT algorithm.
     * 
     * @param dat data to perform FFT on.
     * @return vector<complex<double>> transformed data.
     */
    vector<complex<double>> cooley_tukey(vector<complex<double>> &dat)
    {

        const double pi = acos(-1.0);

        // Get size of input vector.
        uint64_t N = dat.size();
        if (N == 1)
        {
            //do nothing
        }
        else
        {
            // Divide N by 2.
            double N_2 = (double)N / 2;

            // Initialize to store DFTs for even and odd indices.
            vector<complex<double>> even;
            vector<complex<double>> odd;

            // Separate input vector into even and odd indices.
            for (uint64_t i = 0; i < (uint64_t)N_2; i++)
            {
                even.push_back(dat[2 * i]);
                odd.push_back(dat[2 * i + 1]);
            }

            // Recursively compute even and odd vectors.
            cooley_tukey(even);
            cooley_tukey(odd);

            complex<double> q = 0;

            for (uint64_t k = 0; k < (uint64_t)N_2; k++)
            {
                // Twiddle factor times odd values.
                q = polar(1.0, (-2.0 * pi * (double)k) / (double)(N)) * odd[k];

                // Update k-th DFT.
                dat[k] = even[k] + q;
                dat[k + (uint64_t)N_2] = even[k] - q;
            }
        }

        return dat;
    }

    /**
     * @brief Simple low pass filter to remove freqencies above a threshold.
     * 
     * @param fft_values FFT data.
     * @return vector<complex<double>> 
     */
    vector<complex<double>> low_pass(vector<complex<double>> &fft_values, uint32_t &freq, uint32_t &f_samp)
    {
        uint32_t N = (uint32_t)fft_values.size();

        // Compute threshold index
        double threshold = round((N * freq) / f_samp);

        // Remove frequency domain outside of threshold
        for (uint32_t i = (uint32_t)threshold; i < ((uint32_t)fft_values.size() - (uint32_t)threshold); i++)
        {
            fft_values[i] = 0;
        }

        return fft_values;
    }

    /**
     * @brief Get the data object.
     * 
     * @return vector<complex<double>> 
     */
    vector<complex<double>> get_data()
    {
        return x;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief FFT values.
     * 
     */
    vector<complex<double>> x;
};

/**
 * @brief Inverse Fast Fourier transform class.
 * 
 */
class i_fft
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new inverse fft object.
     * 
     * @param dft values to perform inverse FFT on.
     */
    i_fft(vector<complex<double>> &dft)
    {
        y = inverse(dft);
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Inverse Cooley-Tukey radix-2 FFT.
     * 
     * @param freq_domain frequency domain values to perform inverse FFT on.
     * @return vector<complex<double>> inverse FFT values.
     */
    vector<complex<double>> inverse(vector<complex<double>> &freq_domain)
    {
        // Compute complex conjugate of data in freqency domain.
        for (uint64_t i = 0; i < freq_domain.size(); i++)
        {
            freq_domain[i] = conj(freq_domain[i]);
        }

        // Compute forward FFT.
        fft out(freq_domain);
        freq_domain = out.get_data();

        uint64_t N = freq_domain.size();

        // Compute complex conjugate and scale output by N.
        for (uint64_t i = 0; i < freq_domain.size(); i++)
        {
            freq_domain[i] = conj(freq_domain[i]) / (double)N;
        }

        return freq_domain;
    }

    /**
     * @brief Write Inverse FFT values to WAV file.
     * 
     * @param header WAV file header.
     */
    int16_t write_wav(string &header)
    {
        stringstream ss_wav;

        // Create output file name.
        chrono::system_clock now;
        time_t output_time = chrono::system_clock::to_time_t(now.now());
        ss_wav << put_time(localtime(&output_time), "%y%m%d_%H%M%S_.wav");

        // Create output WAV file.
        ofstream output_wav(ss_wav.str(), ios::binary);

        if (!output_wav.is_open())
        {
            cout << "Error opening output file!";
            return -1;
        }

        // Write WAV file header.
        for (uint64_t i = 0; i < header.size(); i++)
        {
            output_wav.write(&header[i], sizeof(int8_t));
        }

        // Write WAV data with real part of inverse FFT values.
        for (uint64_t j = 0; j < y.size(); j++)
        {
            double real_part = real(y[j]);
            int16_t round_real = (int16_t)round(real_part);
            output_wav.write((char *)&round_real, sizeof(int16_t));
        }

        output_wav.close();
        return 0;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Inverse Fast Fourier transform values.
     * 
     */
    vector<complex<double>> y;
};

/**
 * @brief frequency class
 * 
 */
class frequency
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new frequency object.
     * 
     * @param num string to convert to frequency class.
     */
    frequency(const string &num)
    {
        uint32_t input_freq = 0;

        // Verify input is an integer.
        input_freq = stol(num);

        // Verify frequency is in audible range.
        if (input_freq < 20 || input_freq > 20000)
        {
            throw out_of_range("Frequency outside audible range (20 Hz, 20000 Hz) \n");
        }

        f = input_freq;
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Get the data object.
     * 
     * @return uint32_t frequency.
     */
    uint32_t get_data()
    {
        return f;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief frequency in Hz.
     * 
     */
    uint32_t f;
};

int main(int argc, char *argv[])
{
    if (argc == 4)
    {
        try
        {
            vector<complex<double>> input_data;
            vector<complex<double>> output_data;
            vector<complex<double>> freq_data;
            string arg_two;

            // Read and validate WAV file.
            wav_file input(argv[1]);

            input_data = input.get_data();

            // Compute FFT on input WAV file.
            fft dft(input_data);
            freq_data = dft.get_data();

            arg_two = argv[2];
            if (arg_two != "low")
            {
                cout << "Invalid argument for filter. Valid arguments are: \"low\" \n";
                return -1;
            }

            frequency input_freq(argv[3]);
            uint32_t freq = input_freq.get_data();
            uint32_t sampling_rate = input.get_samp();

            // Compute low pass filter on frequency domain
            output_data = dft.low_pass(freq_data, freq, sampling_rate);

            // Compute inverse FFT on input FFT values.
            i_fft output(output_data);

            // Write inverse FFT values to output WAV file.
            string input_header = input.get_header();
            int16_t write_status = output.write_wav(input_header);

            if (write_status < 0)
            {
                return -1;
            }
        }

        catch (const invalid_argument &e)
        {
            cout << "Error: " << e.what() << '\n';
            return -1;
        }
        catch (const out_of_range &e)
        {
            cout << "Error: " << e.what() << '\n';
            return -1;
        }
        catch (const ios_base::failure &e)
        {
            cout << "Error: " << e.what() << '\n';
            return -1;
        }
    }
    else
    {
        cout << (argc - 1) << " argument(s) provided to program. Program requires three arguments: \n";
        cout << "1. The name of a WAV audio file: \"filename.wav\" \n";
        cout << "2. Digital filter to apply: \"low\" \n";
        cout << "3. Frequency (Hz) to apply filter with: ex. 800 \n";
    }
}