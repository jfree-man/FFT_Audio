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
#include <algorithm>
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

        string s;              // String to store input.
        string wav_identifier; // String to store expected WAV identifier.
        string format_size;
        string current;            // Current input value.
        int16_t temp;              // Temporary variable to store input values.
        int16_t bits_persamp;      // String to store bits per sample.
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

        input.seekg(16);

        for (uint32_t i = 0; i < 4; i++)
        {
            current = (char)input.get();
            format_size.insert(i, current);
        }

        // Read in WAV header (44 bytes).
        input.seekg(0);

        for (uint32_t k = 0; k < (20 + *(uint32_t *)&format_size[0] + 8); k++)
        {
            current = (char)input.get();
            header.insert(k, current);
        }

        bits_persamp = *(int16_t *)&header[34];
        if (bits_persamp != 16)
        {
            throw invalid_argument("Not a 16-bit .WAV file.");
        }
        // Get data section size.
        data_section = *(int32_t *)&header[20 + *(uint32_t *)&format_size[0] + 4];

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
    void write_wav(string &outfile)
    {
        stringstream ss_wav;

        // Create output file name.
        // chrono::system_clock now;
        // time_t output_time = chrono::system_clock::to_time_t(now.now());
        // ss_wav << put_time(localtime(&output_time), "%y%m%d_%H%M%S_.wav");

        // Create output WAV file.
        ofstream output_wav(outfile, ios::binary);

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

    uint32_t orig_size()
    {
        uint32_t format_size = *(uint32_t *)&header[16];
        uint32_t data_section = *(int32_t *)&header[20 + format_size + 4];
        return data_section;
    }

    void update_size(uint32_t s)
    {
        uint32_t format_size = *(uint32_t *)&header[16];

        header[20 + format_size + 4] = (unsigned char)s;
        header[20 + format_size + 5] = (unsigned char)(s >> 8);
        header[20 + format_size + 6] = (unsigned char)(s >> 16);
        header[20 + format_size + 7] = (unsigned char)(s >> 24);

        uint32_t file_size = (header.size() - 8) + s;

        header[4] = (unsigned char)file_size;
        header[5] = (unsigned char)(file_size >> 8);
        header[6] = (unsigned char)(file_size >> 16);
        header[7] = (unsigned char)(file_size >> 24);
    }
    void update_fsamp(uint32_t f)
    {
        header[24] = (unsigned char)f;
        header[25] = (unsigned char)(f >> 8);
        header[26] = (unsigned char)(f >> 16);
        header[27] = (unsigned char)(f >> 24);
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
        try
        {
            input_freq = stol(num);
        }
        catch (const invalid_argument &e)
        {
            throw invalid_argument("Invalid frequency argument. Expected a number.");
        }
        catch (const invalid_argument &e)
        {
            throw out_of_range("Frequency input is out of range for 'std::stol'.");
        }

        // Verify frequency is in audible range.
        if (input_freq < 20 || input_freq > 20000)
        {
            throw out_of_range("Frequency outside audible range (20 Hz, 20000 Hz). \n");
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
        uint64_t N = x.size();
        uint32_t freq = f.get_frequency();

        //double i = round((double)freq / ((double)f_samp / (double)N));
        double i = round((double)freq * (double)N) / ((double)f_samp);
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
        threshold_one = index;
        bins = fft_values;
        n = (uint32_t)fft_values.size();
    }
    //constructor
    filter(uint32_t &index_one, uint32_t &index_two, vector<complex<double>> &fft_values)
    {
        //should index, be fft_values.size()/2
        if ((index_one > fft_values.size()) && (index_two > fft_values.size()))
        {
            throw out_of_range("Filter index value is larger than FFT vector");
        }
        threshold_one = index_one;
        threshold_two = index_two;
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
        for (uint32_t i = threshold_one; i < (n - threshold_one); i++)
        {
            bins[i] = 0;
        }
        return bins;
    }
    vector<complex<double>> high_pass()
    {
        // Remove frequency domain outside of threshold
        for (uint32_t i = 0; i < threshold_one; i++)
        {
            bins[i] = 0;
        }
        for (uint32_t i = (n - threshold_one); i < n; i++)
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

    vector<complex<double>> band_pass()
    {
        for (uint32_t i = 0; i < threshold_one; i++)
        {
            bins[i] = 0;
        }
        for (uint32_t i = threshold_two; i < (n / 2 + threshold_two); i++)
        {
            bins[i] = 0;
        }

        for (uint32_t i = (n - threshold_one); i < n; i++)
        {
            bins[i] = 0;
        }

        return bins;
    }

private:
    uint32_t threshold_one;
    uint32_t threshold_two;
    vector<complex<double>> bins;
    uint32_t n;
};
class equalize
{

public:
    //constructor
    equalize()
    {
        modify.assign({0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
        bands = get_bands();
    }
    equalize(char *one, char *two, char *three, char *four, char *five, char *six, char *seven, char *eight, char *nine, char *ten)
    {
        //input decibels
        vector<string> input_vec = {one, two, three, four, five, six, seven, eight, nine, ten};
        int32_t input_dec;
        for (uint32_t i = 0; i < input_vec.size(); i++)
        {
            //verify input is an integer
            try
            {
                input_dec = stol(input_vec[i]);
            }
            catch (const invalid_argument &e)
            {
                throw invalid_argument("Invalid equalizing scale argument. Expected a number. \n");
            }
            catch (const invalid_argument &e)
            {
                throw out_of_range("Equalizing scaler is out of range for 'std::stol'. \n");
            }

            //verify input is in decibel range (-12-12)?
            if (input_dec < -24 || input_dec > 24)
            {
                throw out_of_range("Equalizing scale input is outside valid range. \n");
            }

            modify.push_back(input_dec);
        }

        bands = get_bands();
    }

    //flat band
    vector<complex<double>> flat_band(fft &fft, uint32_t &fs)
    {

        //need freq indices
        //uint32_t lower_index;
        //uint32_t upper_index;
        vector<complex<double>> fft_data;
        vector<uint32_t> indices;
        double m;
        double gain;

        fft_data = fft.get_data();

        //get indices first
        for (uint32_t i = 0; i < bands.size(); i++)
        {
            indices.push_back(fft.get_index(bands[i], fs));
        }

        //for each band
        //if not 0
        for (uint32_t i = 0; i < (bands.size() - 1); i++)
        {
            gain = pow(10.0, ((double)modify[i] / 20.0));

            for (uint32_t k = indices[i]; k < indices[i + 1]; k++)
            {
                fft_data[k] = fft_data[k] * gain;
            }
        }

        return fft_data;
    }

    vector<frequency> get_bands()
    {
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

        bands.assign({f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11});

        return bands;
    }

    //peak band
private:
    vector<int32_t> modify;
    //these should be frequencies
    vector<frequency> bands;
};

class process_audio
{
public:
    process_audio(wav_file &input, string &type, frequency &f1, frequency &f2, equalize &e)
    {
        process(input, type, f1, f2, e);
    }

    process_audio(wav_file &input_one, wav_file &input_two, string &type)
    {
        if (input_one.get_channel() != input_two.get_channel())
        {
            throw invalid_argument("Input files have different number of channels. Cannot mix files.\n");
        }

        if (input_one.get_samp() != input_two.get_samp())
        {
            cout << "Input files have different sampling rates. Using larger sampling rate. \n";

            if (input_two.get_samp() > input_one.get_samp())
            {
                uint32_t f = input_two.get_samp();
                input_one.update_fsamp(f);
            }
        }

        if (type == "add")
        {
            add(input_one, input_two);
        }

        if (type == "overlap")
        {
            overlap(input_one, input_two);
        }
    }

    void add(wav_file &input_one, wav_file &input_two)
    {
        vector<complex<double>> data_one;
        vector<complex<double>> data_two;

        data_one = input_one.get_data();
        data_two = input_two.get_data();

        uint32_t s = input_one.orig_size() + input_two.orig_size();
        uint32_t input_one_bitsize = input_one.orig_size() / 2;
        uint32_t input_two_bitsize = input_two.orig_size() / 2;

        input_one.update_size(s);

        output_data = data_one;

        //concatenate files
        output_data.insert(output_data.begin() + input_one_bitsize, data_two.begin(), data_two.begin() + input_two_bitsize);

        //update header with sampling rate
    }

    void overlap(wav_file &input_one, wav_file &input_two)
    {
        vector<complex<double>> data_one;
        vector<complex<double>> data_two;
        uint32_t m; //max size
        uint32_t maximum_flag = 0;
        complex<double> max_value = 0;

        data_one = input_one.get_data();
        data_two = input_two.get_data();

        uint32_t input_one_bitsize = input_one.orig_size() / 2;
        uint32_t input_two_bitsize = input_two.orig_size() / 2;

        if (input_one_bitsize > input_two_bitsize)
        {
            for (uint32_t i = 0; i < input_one_bitsize; i++)
            {
                if (i < input_two_bitsize)
                {
                    output_data.push_back(data_one[i] + data_two[i]);
                }
                else
                {
                    output_data.push_back(data_one[i]);
                }

                if (real(output_data[i]) > real(max_value))
                {
                    max_value = output_data[i];
                }
            }
        }
        else
        {
            //input two is larger
            m = input_two.orig_size();
            input_one.update_size(m);
            for (uint32_t i = 0; i < input_two_bitsize; i++)
            {
                if (i < input_one_bitsize)
                {
                    output_data.push_back(data_one[i] + data_two[i]);
                }
                else
                {
                    output_data.push_back(data_two[i]);
                }

                if (real(output_data[i]) > real(max_value))
                {
                    max_value = output_data[i];
                }
            }
        }

        if (real(max_value) > 32768)
        {
            double scale = real(max_value) / 32768;
            for (uint32_t i = 0; i < output_data.size(); i++)
            {
                output_data[i] = round(real(output_data[i]) / scale) + 0i;
            }
        }
    }

    void process(wav_file &input, string &type, frequency &f1, frequency &f2, equalize &e)
    {
        cout << "Processing WAV file...\n";
        vector<complex<double>> input_data;

        vector<complex<double>> fft_data;
        vector<complex<double>> ifft_data;
        vector<complex<double>> freq_data;
        vector<complex<double>> transformed_data;
        uint32_t sampling_rate;

        input_data = input.get_data();
        uint64_t input_size = input_data.size();

        uint16_t channel = input.get_channel();

        sampling_rate = input.get_samp();

        if (channel == 1) //mono
        {
            if (input_size <= 262144)
            {
                //doo processing all in window
                //fft
                fft fft_all(input_data);
                fft_data = fft_all.get_data();

                //process

                //uint32_t index = 1000;
                //100 cut out noises for low pass
                //110 let noises for high pass, 150000 cuts it out, 50000 small noises, 30000 somewhat recognizable,
                //10000 sounds about right, 5000 better
                //filter fil(index, fft_data);
                //transformed_data = fil.low_pass();
                transformed_data = process_type(fft_all, type, f1, f2, e, sampling_rate);

                //ifft
                i_fft ifft_all(transformed_data);

                //output
                output_data = ifft_all.get_data();
            }
            else
            {
                //process in separate windows
                output_data = window(input_size, input_data, type, f1, f2, e, sampling_rate);
            }
        }
        else //channel ==2
        {
            //split up input_data into two left and right channels
            //if each is different

            vector<complex<double>> left;
            vector<complex<double>> right;
            vector<complex<double>> fft_left_data;
            vector<complex<double>> fft_right_data;
            vector<complex<double>> transformed_left;
            vector<complex<double>> transformed_right;
            vector<complex<double>> output_left;
            vector<complex<double>> output_right;
            uint64_t channel_size = input_size / 2;

            for (uint32_t j = 0; j < channel_size; j++)
            {
                left.push_back(input_data[j * 2]);
                right.push_back(input_data[(j * 2) + 1]);
            }
            if (channel_size <= 262144)
            {
                //doo processing all in window
                //fft
                fft fft_left(left);
                fft fft_right(right);
                fft_left_data = fft_left.get_data();
                fft_right_data = fft_right.get_data();

                //process

                //uint32_t index = 8000;
                //100 cut out noises for low pass
                //110 let noises for high pass, 150000 cuts it out, 50000 small noises, 30000 somewhat recognizable,
                // //10000 sounds about right, 5000 better
                // filter fil_left(index, fft_left_data);
                // filter fil_right(index, fft_right_data);
                // transformed_left = fil_left.low_pass();
                // transformed_right = fil_right.low_pass();

                transformed_left = process_type(fft_left, type, f1, f2, e, sampling_rate);
                transformed_right = process_type(fft_right, type, f1, f2, e, sampling_rate);

                //ifft
                i_fft ifft_left(transformed_left);
                i_fft ifft_right(transformed_right);

                //output
                output_left = ifft_left.get_data();
                output_right = ifft_right.get_data();
            }
            else
            {

                //have to window
                //process in separate windows
                output_left = window(channel_size, left, type, f1, f2, e, sampling_rate);
                output_right = window(channel_size, right, type, f1, f2, e, sampling_rate);
            }
            //combine output
            for (uint32_t j = 0; j < channel_size; j++)
            {
                output_data.push_back(output_left[j]);
                output_data.push_back(output_right[j]);
            }
        }
    }

    //member functions

    vector<complex<double>> window(uint64_t &input_size, vector<complex<double>> &input_data, string &type, frequency &f1, frequency &f2, equalize &e, uint32_t &sampling_rate)
    {
        vector<double> ramp; //linear ramp from 0 to 1 of length of window size
        vector<complex<double>> fft_data;
        vector<complex<double>> ifft_data;
        vector<complex<double>> transformed_data;
        vector<complex<double>> final;
        // Window
        //uint64_t window_size = 1024;
        //uint64_t window_size = 2097152;
        uint64_t window_size = 262144;
        uint64_t num_windows = input_size / window_size;
        uint32_t progress_counter = 1; //progress counter for windowing

        ramp.push_back(0);

        //create ramp
        for (uint32_t j = 1; j < (window_size - 1); j++)
        {
            ramp.push_back(((double)j) / ((double)window_size - 1));
        }

        //end ramp
        ramp.push_back(1);

        for (uint32_t i = 0; i < (num_windows - 1); i++)
        {

            cout << "[";
            //number of loops completed
            for (uint32_t progress = 0; progress < progress_counter; progress++)
            {
                cout << "-";
            }
            //number of loops remaining
            for (uint32_t progress = progress_counter; progress < (num_windows - 1); progress++)
            {
                cout << " ";
            }
            //percentage of loop completion
            cout << "]" << round((double)progress_counter / (double)(num_windows - 1) * (double)100) << "% \r";

            //bad to initialize each time?
            vector<complex<double>> chunk_one(input_data.begin() + (i * window_size), input_data.begin() + ((i * window_size) + window_size));
            vector<complex<double>> chunk_two(input_data.begin() + ((i * window_size) + window_size), input_data.begin() + ((i * window_size) + (2 * window_size)));

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

            //frequency input_freq(argv[3]);
            // uint32_t freq = input_freq.get_data();
            //uint32_t sampling_rate = input.get_samp();

            //output_data1 = df1.low_pass(df1_data, freq, sampling_rate);
            // // output_data2 = df2.low_pass(df2_data, freq, sampling_rate);
            ////////////////////////////////////////////////////////////
            //low pass filter
            //uint32_t index = fft_chunk_one.get_index(input_freq, sampling_rate);
            //uint32_t index = 50000;
            //100 cut out noises for low pass
            //110 let noises for high pass, 150000 cuts it out, 50000 small noises, 30000 somewhat recognizable,
            //10000 sounds about right, 5000 better
            //filter fil(index, fft_data);
            //transformed_data = fil.low_pass();
            //transformed_data = fil.high_pass();
            transformed_data = process_type(fft_chunk_one, type, f1, f2, e, sampling_rate);
            ///////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////
            //testing
            //i_fft if1(output_data1);
            i_fft ifft_chunk_one(transformed_data);

            //testing

            ifft_data = ifft_chunk_one.get_data();

            if (i == 0)
            {
                final.insert(final.begin(), input_data.size(), 0);
            }

            for (uint32_t k = 0; k < ifft_data.size(); k++)
            {
                //starting position

                final[(i * window_size) + k] = final[(i * window_size) + k] + ifft_data[k];
                //cout << "output" << output[(i * window_size) + k] << "\n";
            }

            progress_counter++;
        }
        cout << "\n";
        return final;
    }

    vector<complex<double>> process_type(fft &_fft, string &type, frequency &f1, frequency &f2, equalize &e, uint32_t f_samp)
    {
        uint32_t index;
        uint32_t index_two;
        vector<complex<double>> fft_data;
        vector<complex<double>> transformed_data;
        fft_data = _fft.get_data();
        if (type == "low")
        {
            //get index

            index = _fft.get_index(f1, f_samp);
            //declare filter
            filter low_fil(index, fft_data);

            //get member function
            transformed_data = low_fil.low_pass();
        }
        if (type == "high")
        {
            //get index
            index = _fft.get_index(f1, f_samp);

            //declare filter
            filter high_fil(index, fft_data);

            //get member function
            transformed_data = high_fil.high_pass();
        }
        if (type == "band")
        {
            //add this
            index = _fft.get_index(f1, f_samp);
            index_two = _fft.get_index(f2, f_samp);

            filter band_fil(index, index_two, fft_data);

            transformed_data = band_fil.band_pass();
        }
        if (type == "equalize")
        {
            transformed_data = e.flat_band(_fft, f_samp);
        }
        return transformed_data;
    }

    vector<complex<double>> get_data()
    {
        return output_data;
    }

private:
    vector<complex<double>> output_data;
};
int main(int argc, char *argv[])
{
    if (argc == 5 || argc == 6 || argc == 14)
    {
        //for band pass first argument must be smaller

        try
        {
            vector<complex<double>> input_data;
            vector<complex<double>> output_data;
            string arg_two;   //do I need to check outgile name, are there restructions or will any string
            string arg_three; //low, high, band, equalize

            arg_two = argv[2];
            if (arg_two.size() < 4)
            {
                throw invalid_argument("Output file name is less that 4 characters.");
            }
            if (arg_two.substr((arg_two.size() - 4), 4) != ".wav")
            {
                throw invalid_argument("Output file name is missing \".wav\" extension.");
            }

            // Read and validate WAV file.
            wav_file input(argv[1]);
            arg_three = argv[3];

            switch (argc)
            {
            case 5:

                if (arg_three != "low" && arg_three != "high" && arg_three != "add" && arg_three != "overlap")
                {
                    throw invalid_argument("Provided 4 arguments to program. Third argument is not \"low\" or \"high\" or \"add\". \n");
                }
                else if (arg_three == "low" || arg_three == "high")
                {
                    frequency input_freq(argv[4]);
                    equalize buffer;
                    process_audio processor(input, arg_three, input_freq, input_freq, buffer);
                    output_data = processor.get_data();
                }
                else
                {
                    wav_file input_two(argv[4]);
                    process_audio mixer(input, input_two, arg_three);
                    output_data = mixer.get_data();
                }
                break;

            case 6:
                if (arg_three != "band")
                {
                    throw invalid_argument("Provided 5 arguments to program. Third argument is not \"band\". \n");
                }
                else
                {
                    frequency lower_freq(argv[4]);
                    frequency upper_freq(argv[5]);

                    uint32_t lower = lower_freq.get_frequency();
                    uint32_t upper = upper_freq.get_frequency();

                    if (upper <= lower)
                    {
                        throw invalid_argument("Second frequency argument is smaller than first frequency argument. \n");
                    }
                    equalize buffer;
                    process_audio processor(input, arg_three, lower_freq, upper_freq, buffer);
                    output_data = processor.get_data();
                }
                break;
            case 14:
                if (arg_three != "equalize")
                {
                    throw invalid_argument("Provided 13 arguments to program. Third argument is not \"equalize\". \n");
                }
                else
                {
                    uint32_t buffer = 20;
                    frequency buffer_1(buffer);
                    frequency buffer_2(buffer);
                    equalize modify(argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12], argv[13]);
                    process_audio processor(input, arg_three, buffer_1, buffer_2, modify);
                    output_data = processor.get_data();
                }

                break;
            }

            //////////////////FILE WRITING////////////////////////////////////
            // cout << "output data " << output_data.size();
            wav_file output(input, output_data);
            output.write_wav(arg_two);
            //cout << "end\n";
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
        cout << (argc - 1) << " argument(s) provided to program. Program expects: \n";
        cout << "1. The name of a WAV audio file. Example: \"infile.wav\" \n";
        cout << "2. The name of the desired output WAV audio file. Example: \"outfile.wav\" \n";
        cout << "3. The type of audio processing to perform. Valid arguments are: \"low\", \"high\", \"band\", \"equalize\" \n";
        cout << "4. The arguments to be passed to the type of audio process. \n";
        cout << "\t low: One frequency value in (20,20000) Hz to apply the low-pass filter. Example: 1000 \n";
        cout << "\t high: One frequency value in (20,20000) Hz to apply the high-pass filter. Example: 1000 \n";
        cout << "\t band: Two frequency values in (20,20000) Hz to apply the band-pass filter. Example: 1000 6000\n";
        cout << "\t equalize: Ten integer values between (-24,24) to scale frequency bands. Example: 1 12 0 0 -3 8 -8 9 1 0  \n";
    }
}