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

    wav_file(wav_file &wav, vector<complex<double>> &dat)
    {
        //use the same header and sampling info
        header = wav.header;
        f_samp = wav.f_samp;

        //convert vector to int16_t
        for (uint64_t j = 0; j < dat.size(); j++)
        {
            double real_part = real(dat[j]);
            int16_t round_real = (int16_t)round(real_part);
            data.push_back(round_real);
        }
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

        string s;                  // String to store input.
        string wav_identifier;     // String to store expected WAV identifier.
        string current;            // Current input value.
        int16_t temp;              // Temporary variable to store input values.
        uint64_t data_section = 0; // Size of data section in file.
        uint64_t power = 1;        // Nearest power of 2 to data size.
        uint64_t data_size = 1;    // Size of data in 2-bytes.
        //uint64_t max_data_size = 262144; // Maximum allowable data size (2^18).

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

        for (uint32_t k = 0; k < 44; k++)
        {
            current = (char)input.get();
            header.insert(k, current);
        }

        // Get data section size.
        data_section = *(int32_t *)&header[40];

        // Ensure data section is a multiple of 2.
        if (data_section % 2 != 0)
        {
            data_section = data_section + 1;
        }

        // Calculate size of data in 16 bits.
        data_size = data_section / 2;

        // if (data_size > max_data_size)
        // {
        //     throw invalid_argument("File size is larger than maximum allowable.");
        // }

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
        f_samp = *(uint32_t *)&header[24];

        //Get the number of channels
        channel = *(uint16_t *)&header[22];
    }
    /**
     * @brief Write Inverse FFT values to WAV file.
     * 
     * @param header WAV file header.
     */
    void write_wav()
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
            throw ios_base::failure("Error opening input file! \n");
        }

        // Write WAV file header.
        for (uint64_t i = 0; i < header.size(); i++)
        {
            output_wav.write(&header[i], sizeof(int8_t));
        }

        // Write WAV data with real part of inverse FFT values.
        for (uint64_t j = 0; j < data.size(); j++)
        {
            output_wav.write((char *)&data[j], sizeof(int16_t));
        }

        output_wav.close();
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

    uint16_t get_channel()
    {
        return channel;
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

    uint16_t channel;
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

    frequency(uint32_t &num)
    {
        // Verify frequency is in audible range.
        if (num < 20 || num > 20000)
        {
            throw out_of_range("Frequency outside audible range (20 Hz, 20000 Hz) \n");
        }
        f = num;
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Get the data object.
     * 
     * @return uint32_t frequency.
     */
    uint32_t get_frequency()
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
     * @brief Get the data object.
     * 
     * @return vector<complex<double>> 
     */
    vector<complex<double>> get_data()
    {
        return x;
    }

    //should sampling frequency be passed somehow
    uint32_t get_index(frequency &f, uint32_t f_samp)
    {
        uint32_t N = x.size();
        uint32_t freq = f.get_frequency();

        //double i = round((double)freq / ((double)f_samp / (double)N));
        double i = round((double)freq / ((double)f_samp / (double)N));
        uint32_t index = (uint32_t)i;
        return index;
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

    vector<complex<double>> get_data()
    {
        return y;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Inverse Fast Fourier transform values.
     * 
     */
    vector<complex<double>>
        y;
};

class filter
{
public:
    //constructor
    filter(uint32_t &index, vector<complex<double>> &fft_values)
    {
        //should index, be fft_values.size()/2
        if (index > fft_values.size())
        {
            throw out_of_range("Filter index value is larger than FFT vector");
        }
        threshold = index;
        bins = fft_values;
        n = (uint32_t)fft_values.size();
    }

    /**
     * @brief Low-pass sharp-cut digital filter.
     * 
     * @param fft_values vector of FFT values.
     * @param freq frequency to apply filter with.
     * @param f_samp sampling rate.
     * @return vector<complex<double>> filtered FFT values.
     */
    vector<complex<double>> low_pass()
    {
        // Remove frequency domain outside of threshold
        for (uint32_t i = threshold; i < (n - threshold); i++)
        {
            bins[i] = 0;
        }
        return bins;
    }
    vector<complex<double>> high_pass()
    {
        // Remove frequency domain outside of threshold
        for (uint32_t i = 0; i < threshold; i++)
        {
            bins[i] = 0;
        }
        for (uint32_t i = (n - threshold); i < n; i++)
        {
            bins[i] = 0;
        }

        // for (uint32_t i = 0; i < bins.size(); i++)
        // {
        //     cout << bins[i] << "\n";
        // }
        return bins;
    }

    //bandpass

    //notch

private:
    uint32_t threshold;
    vector<complex<double>> bins;
    uint32_t n;
};
class equalize
{

public:
    //constructor
    equalize(string &one, string &two, string &three, string &four, string &five, string &six, string &seven, string &eight, string &nine, string &ten)
    {
        //input decibels
        vector<string> input_vec = {one, two, three, four, five, six, seven, eight, nine, ten};
        uint32_t input_dec;
        for (uint32_t i = 0; i < input_vec.size(); i++)
        {
            //verify input is an integer
            input_dec = stol(input_vec[i]);

            //verify input is in decibel range (-12-12)?
            if (input_dec < -12 || input_dec > 12)
            {
                throw out_of_range("Decibel input is outside valid range. \n");
            }

            modify.push_back(input_dec);
        }

        uint32_t n1 = 23;
        uint32_t n2 = 45;
        uint32_t n3 = 89;
        uint32_t n4 = 177;
        uint32_t n5 = 354;
        uint32_t n6 = 707;
        uint32_t n7 = 1414;
        uint32_t n8 = 2828;
        uint32_t n9 = 5657;
        uint32_t n10 = 11314;
        uint32_t n11 = 20000;

        frequency f1(n1);
        frequency f2(n2);
        frequency f3(n3);
        frequency f4(n4);
        frequency f5(n5);
        frequency f6(n6);
        frequency f7(n7);
        frequency f8(n8);
        frequency f9(n9);
        frequency f10(n10);
        frequency f11(n11);

        bands[0] = f1;
        bands[1] = f2;
        bands[2] = f3;
        bands[3] = f4;
        bands[4] = f5;
        bands[5] = f6;
        bands[6] = f7;
        bands[7] = f8;
        bands[8] = f9;
        bands[9] = f10;
        bands[10] = f11;
    }

    //flat band
    void flat_band(fft &fft, wav_file &wav)
    {
        //get sampling rate
        uint32_t fs = wav.get_samp();

        //need freq indices
        uint32_t lower_index;
        uint32_t upper_index;
        vector<complex<double>> fft_data;
        vector<uint32_t> indices;

        fft_data = fft.get_data();

        //get indices first
        for (uint32_t i = 0; i < (bands.size() - 1); i++)
        {
            indices.push_back(fft.get_index(bands[i], fs));
        }

        //for each band
        //if not 0
        for (uint32_t i = 0; i < (bands.size() - 1); i++)
        {
            if (modify[i] != 0)
            {
                for (uint32_t k = indices[i]; k < indices[i + 1]; k++)
                {
                    fft_data[k] = fft_data[k] * (double)modify[i];
                }
            }
        }
    }

    //peak band
private:
    vector<int32_t> modify;
    //these should be frequencies
    vector<frequency> bands;
};

int main(int argc, char *argv[])
{
    if (argc == 4)
    {
        try
        {
            vector<complex<double>> input_data;
            vector<complex<double>> output_data;
            vector<complex<double>> fft_data;
            vector<complex<double>> ifft_data;
            vector<complex<double>> freq_data;
            vector<complex<double>> transformed_data;
            string arg_two;

            uint16_t channel;

            // Read and validate WAV file.
            wav_file input(argv[1]);

            //Get number of channels.
            channel = input.get_channel();

            input_data = input.get_data();

            // Window
            //uint64_t window_size = 1024;
            //uint64_t window_size = 2097152;
            uint64_t window_size = 262144;
            uint64_t num_windows = input_data.size() / window_size;
            cout << input_data.size() << "\n";
            cout << num_windows << "\n";

            vector<double> ramp; //linear ramp from 0 to 1 of length of window size
            ramp.push_back(0);

            //create ramp
            for (uint32_t j = 1; j < (window_size - 1); j++)
            {
                ramp.push_back(((double)j) / ((double)window_size - 1));
            }

            //end ramp
            ramp.push_back(1);

            //seems right
            // for (uint32_t j = 0; j < ramp.size(); j++)
            // {
            //     cout << "at " << j << ramp[j] << "\n";
            // }

            // separate input data into 2^something parts

            if (channel == 2)
            {
                //do i make below a class and function of its own
            }

            for (uint32_t i = 0; i < (num_windows - 1); i++)
            {
                //bad to initialize each time?
                vector<complex<double>> chunk_one(input_data.begin() + (i * window_size), input_data.begin() + ((i * window_size) + window_size));
                vector<complex<double>> chunk_two(input_data.begin() + ((i * window_size) + window_size), input_data.begin() + ((i * window_size) + (2 * window_size)));
                // for (uint32_t k = 0; k < chunk_one.size(); k++)
                // {
                //     cout << "chunk one" << chunk_one[k] << "\n";
                // }
                // for (uint32_t k = 0; k < chunk_two.size(); k++)
                // {
                //     cout << "chunk two" << chunk_two[k] << "\n";
                // }

                //transform data by linear ramp
                //multiply chunk one by ramp
                //chunk two = chunk two - (ramp*chunk_two)
                ///////////////////////////////////////////////////////
                for (uint32_t k = 0; k < chunk_one.size(); k++)
                {

                    chunk_one[k] = chunk_one[k] * ramp[k];
                    chunk_two[k] = chunk_two[k] - (chunk_two[k] * ramp[k]);
                }

                /////////////////////////////////////////////////////////////

                //fft big chunk
                //make chucnk one big vector
                chunk_one.insert(chunk_one.end(), chunk_two.begin(), chunk_two.end());
                fft fft_chunk_one(chunk_one);

                //do the ifft on each part
                fft_data = fft_chunk_one.get_data();

                frequency input_freq(argv[3]);
                // uint32_t freq = input_freq.get_data();
                uint32_t sampling_rate = input.get_samp();

                //output_data1 = df1.low_pass(df1_data, freq, sampling_rate);
                // // output_data2 = df2.low_pass(df2_data, freq, sampling_rate);
                ////////////////////////////////////////////////////////////
                //low pass filter
                //uint32_t index = fft_chunk_one.get_index(input_freq, sampling_rate);
                uint32_t index = 4800;
                //100 cut out noises for low pass
                //110 let noises for high pass, 150000 cuts it out, 50000 small noises, 30000 somewhat recognizable,
                //10000 sounds about right, 5000 better
                filter fil(index, fft_data);
                transformed_data = fil.low_pass();
                //transformed_data = fil.high_pass();
                ///////////////////////////////////////////////////////

                ////////////////////////////////////////////////////////
                //testing
                //i_fft if1(output_data1);
                i_fft ifft_chunk_one(transformed_data);

                //testing

                ifft_data = ifft_chunk_one.get_data();

                ///////////////////////////////////////////////////////////////

                //output.insert(output.begin(), chunk_one.size(), 0);

                //cout << output.size() << "\n";

                // for (uint32_t k = 0; k < chunk_one.size(); k++)
                // {
                //     output[k] = dat1[k] + dat2[k];
                // }
                // cout << output.size() << "\n";

                //add ifft in specific way
                // ///////////////////////////////////////////////////////////////
                if (i == 0)
                {
                    output_data.insert(output_data.begin(), input_data.size(), 0);
                }

                /////////////////
                // for (uint32_t k = 0; k < chunk_one.size(); k++)
                // {
                //     //starting position
                //     output[(i * window_size) + k] = output[(i * window_size) + k] + dat1[k];
                // }
                // for (uint32_t k = 0; k < chunk_two.size(); k++)
                // {
                //     //starting position
                //     output[(i * window_size) + window_size + k] = output[(i * window_size) + window_size + k] + dat2[k];
                // }
                // ///////////////////////////////////////////////////////////////

                for (uint32_t k = 0; k < ifft_data.size(); k++)
                {
                    //starting position

                    output_data[(i * window_size) + k] = output_data[(i * window_size) + k] + ifft_data[k];
                    //cout << "output" << output[(i * window_size) + k] << "\n";
                }

                cout << "loop end \n";
            }

            // for (uint32_t k = 0; k < output.size(); k++)
            // {
            //     cout << "final output: " << output[k] << "\n";
            // }
            // // Compute FFT on input WAV file.
            // fft dft(input_data);
            // freq_data = dft.get_data();

            // arg_two = argv[2];
            // if (arg_two != "low")
            // {
            //     cout << "Invalid argument for filter. Valid arguments are: \"low\" \n";
            //     return -1;
            // }

            // frequency input_freq(argv[3]);
            // uint32_t freq = input_freq.get_data();
            // uint32_t sampling_rate = input.get_samp();

            // // Compute low pass filter on frequency domain
            // output_data = dft.low_pass(freq_data, freq, sampling_rate);

            // // Compute inverse FFT on input FFT values.
            // i_fft output(output_data);

            // Write inverse FFT values to output WAV file.
            //string input_header = input.get_header();
            //int16_t write_status = output.write_wav(input_header);

            // if (write_status < 0)
            // {
            //     return -1;
            // }

            //////////////////TEMP FILE WRITING////////////////////////////////////
            wav_file output(input, output_data);
            output.write_wav();
            cout << "end\n";
            //////////////////////////////////////////////////////////////////
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