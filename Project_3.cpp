/**
 * @file main.cpp
 * @author Jennifer Freeman (freemjc@mcmaster.ca)
 * @brief 
 * @version 0.1
 * @date 2021-12-19
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

// ============================================================================================= //
//                                    Begin class wav_file                                       //

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
      * @brief Construct a new wav file object from file name.
      * 
      * @param input_string string of wav file name, (ex. "filename.wav").
      */
    wav_file(const string &input_string)
    {
        read_wav(input_string);
    }

    /**
     * @brief Construct a new wav file object from existing wav_file object.
     * 
     * @param wav wav_file object to extract header and sampling rate.
     * @param dat vector containing audio signal data.
     */
    wav_file(wav_file &wav, vector<complex<double>> &dat)
    {
        // Use the same header and sampling rate from input wav_file.
        header = wav.header;
        f_samp = wav.f_samp;

        // Convert data to 16-bit.
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

        string s;                       // String to store input.
        string wav_identifier;          // String to store expected WAV identifier.
        string format_size;             // Size of format in header.
        string current;                 // Current input value.
        int16_t temp;                   // Temporary variable to store input values.
        int16_t bits_persamp;           // Bits per sample.
        uint64_t data_section = 0;      // Size of data section in file.
        uint64_t power = 1;             // Nearest power of 2 to data size.
        uint64_t data_size = 1;         // Size of data in 2-bytes.
        uint64_t max_data_size = 1.2e7; // Maximum allowable data size (12 MB).

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
        // Get size of format data in header.
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

        // Verify file is 16-bit.
        bits_persamp = *(int16_t *)&header[34];
        if (bits_persamp != 16)
        {
            throw invalid_argument("Not a 16-bit .WAV file.");
        }

        // Get data section size.
        data_section = *(int32_t *)&header[20 + *(uint32_t *)&format_size[0] + 4];

        // Check file size.
        if (data_section > max_data_size)
        {
            throw invalid_argument("File size is larger than maximum allowable (12 MB).");
        }

        // Ensure data section is a multiple of 2.
        if (data_section % 2 != 0)
        {
            data_section = data_section + 1;
        }

        // Calculate size of data in 16 bits.
        data_size = data_section / 2;

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
     * @brief Write WAV file.
     * 
     * @param outfile string of outfile name (ex. "Outfile.wav")
     */
    void write_wav(string &outfile)
    {
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

        // Write WAV data with real part of complex inverse FFT values.
        for (uint64_t j = 0; j < data.size(); j++)
        {
            output_wav.write((char *)&data[j], sizeof(int16_t));
        }

        output_wav.close();
    }

    /**
     * @brief Size of padded data in WAV file.
     * 
     * @return uint64_t size of vector.
     */
    uint64_t datasize()
    {
        return data.size();
    }

    /**
     * @brief Get the 16-bit audio signal data and convert to complex<double>.
     * 
     * @return vector<complex<double>> vector of converted signal.
     */
    vector<complex<double>> get_data()
    {
        vector<complex<double>> new_vec(data.begin(), data.end());
        return new_vec;
    }

    /**
     * @brief Get the 44-byte header string from the WAV file.
     * 
     * @return string header string.
     */
    string get_header()
    {
        return header;
    }

    /**
     * @brief Get the sampling rate from the WAV file..
     * 
     * @return uint32_t sampling rate (ex. 44100 Hz)
     */
    uint32_t get_samp()
    {
        return f_samp;
    }

    /**
     * @brief Get the channel number of the WAV file.
     * 
     * @return uint16_t channel number. 1 = mono, 2 = stereo.
     */
    uint16_t get_channel()
    {
        return channel;
    }

    /**
     * @brief Get the size of the data section from the WAV file in bytes (excluding padding).
     * 
     * @return uint32_t size of the data section in bytes.
     */
    uint32_t orig_size()
    {
        uint32_t format_size = *(uint32_t *)&header[16];
        uint32_t data_section = *(int32_t *)&header[20 + format_size + 4];
        return data_section;
    }

    /**
     * @brief Update header to reflect new WAV file sizes.
     * 
     * @param s size in bytes of data section in new WAV file.
     */
    void update_size(uint32_t s)
    {
        uint32_t format_size = *(uint32_t *)&header[16];

        // Update data section size.
        header[20 + format_size + 4] = (unsigned char)s;
        header[20 + format_size + 5] = (unsigned char)(s >> 8);
        header[20 + format_size + 6] = (unsigned char)(s >> 16);
        header[20 + format_size + 7] = (unsigned char)(s >> 24);

        uint32_t file_size = ((uint32_t)header.size() - 8) + s;

        // Update total file size.
        header[4] = (unsigned char)file_size;
        header[5] = (unsigned char)(file_size >> 8);
        header[6] = (unsigned char)(file_size >> 16);
        header[7] = (unsigned char)(file_size >> 24);
    }

    /**
     * @brief Update header to reflect new sampling rate for new WAV files.
     * 
     * @param f sampling rate.
     */
    void update_fsamp(uint32_t f)
    {
        // Update sampling rate in header.
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
     * @brief Sampling rate in Hz.
     * 
     */
    uint32_t f_samp;

    /**
     * @brief Number of channels in WAV file.
     * 
     */
    uint16_t channel;
};

//                                     End class wav_file                                        //
// ============================================================================================= //

// ============================================================================================= //
//                                   Begin class frequency                                       //

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
     * @brief Construct frequency value from a string.
     * 
     * @param num string to convert to frequency.
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
        catch (const out_of_range &e)
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

    /**
     * @brief Construct a frequency value from a number.
     * 
     * @param num number to convert to frequency.
     */
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
     * @brief Get the frequency value.
     * 
     * @return uint32_t frequency (Hz).
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
//                                     End class frequency                                       //
// ============================================================================================= //
// ============================================================================================= //
//                                       Begin class fft                                         //

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
     * @param dat vector of complex data to perform FFT on.
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
     * @return vector<complex<double>> FFT bins.
     */
    vector<complex<double>> cooley_tukey(vector<complex<double>> &dat)
    {

        const double pi = acos(-1.0);

        // Get size of input vector.
        uint64_t N = dat.size();
        if (N == 1)
        {
            // Do nothing.
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
     * @brief Get the output from the FFT.
     * 
     * @return vector<complex<double>> FFT bins.
     */
    vector<complex<double>> get_data()
    {
        return x;
    }

    /**
     * @brief Get the index of the FFT bins at a given frequency.
     * 
     * @param f frequency.
     * @param f_samp sampling rate of signal.
     * @return uint32_t index of FFT bin.
     */
    uint32_t get_index(frequency &f, uint32_t f_samp)
    {
        // Size of FFT chunk.
        uint64_t N = x.size();

        // Frequency value.
        uint32_t freq = f.get_frequency();

        // Comput index of FFT bin at specified frequency.
        double i = round((double)freq * (double)N) / ((double)f_samp);
        uint32_t index = (uint32_t)i;
        return index;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief values from output of FFT (FFT bins).
     * 
     */
    vector<complex<double>> x;
};

//                                     End class fft                                             //
// ============================================================================================= //

// ============================================================================================= //
//                                   Begin class i_fft                                           //

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
     * @brief Get the inverse FFT values.
     * 
     * @return vector<complex<double>> values from the output of the inverse FFT.
     */
    vector<complex<double>> get_data()
    {
        return y;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Values from the output of the inverse FFT.
     * 
     */
    vector<complex<double>>
        y;
};
//                                     End class i_fft                                           //
// ============================================================================================= //

// ============================================================================================= //
//                                    Begin class filter                                         //
/**
 * @brief Digital filter class for signal processing.
 * 
 */
class filter
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new filter object with one index to apply to FFT bins.
     * 
     * @param index index of frequency value used in filter.
     * @param fft_values FFT bins to filter.
     */
    filter(uint32_t &index, vector<complex<double>> &fft_values)
    {
        // Check index is not larger than length of FFT bins.
        if (index > fft_values.size())
        {
            throw out_of_range("Filter index value is larger than FFT vector");
        }
        threshold_one = index;
        bins = fft_values;
        n = (uint32_t)fft_values.size();
    }

    /**
     * @brief Construct a new filter object with two indicdes to apply to FFT bins.
     * 
     * @param index_one index of frequency value used in filter.
     * @param index_two second index of frequency value used in filter.
     * @param fft_values FFT bins to filter.
     */
    filter(uint32_t &index_one, uint32_t &index_two, vector<complex<double>> &fft_values)
    {
        // Check index is not larger than length of FFT bins.
        if ((index_one > fft_values.size()) && (index_two > fft_values.size()))
        {
            throw out_of_range("Filter index value is larger than FFT vector");
        }
        threshold_one = index_one;
        threshold_two = index_two;
        bins = fft_values;
        n = (uint32_t)fft_values.size();
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Low-pass sharp-cut digital filter. Frequencies above threshold will be removed.
     * 
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

    /**
     * @brief High-pass sharp-cut digital filter. Frequencies below threshold will be removed.
     * 
     * @return vector<complex<double>> filtered FFT values.
     */
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

        return bins;
    }

    /**
     * @brief Band-pass sharp-cut digital filter. Frequencies above and below both thresholds will be removed.
     * 
     * @return vector<complex<double>> filtered FFT values.
     */
    vector<complex<double>> band_pass()
    {
        // Remove frequency domain outside of both thresholds.
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
    // ============
    // Private data
    // ============

    /**
     * @brief Index of first frequency value in filter.
     * 
     */
    uint32_t threshold_one;

    /**
     * @brief Index of second frequency value in filter.
     * 
     */
    uint32_t threshold_two;

    /**
     * @brief Filtered FFT bins.
     * 
     */
    vector<complex<double>> bins;

    /**
     * @brief Length of FFT bins.
     * 
     */
    uint32_t n;
};

//                                     End class filter                                          //
// ============================================================================================= //

// ============================================================================================= //
//                                    Begin class equalize                                       //
/**
 * @brief Equalizer class to perform decibel attenuation and boost to signal.
 * 
 */
class equalize
{

public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new equalize object with default values.
     * 
     */
    equalize()
    {
        modify.assign({0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
        bands = get_bands();
    }

    /**
     * @brief Construct a new equalize object given ten input decibel changes to apply to the ten frequency bands.
     * 
     * @param one decibel change to apply to band one.
     * @param two decibel change to apply to band two.
     * @param three decibel change to apply to band three.
     * @param four decibel change to apply to band four.
     * @param five decibel change to apply to band five.
     * @param six decibel change to apply to band six.
     * @param seven decibel change to apply to band seven.
     * @param eight decibel change to apply to band eight.
     * @param nine decibel change to apply to band nine.
     * @param ten decibel change to apply to band ten.
     */
    equalize(char *one, char *two, char *three, char *four, char *five, char *six, char *seven, char *eight, char *nine, char *ten)
    {
        // Put input decibel changes in a vector.
        vector<string> input_vec = {one, two, three, four, five, six, seven, eight, nine, ten};
        int32_t input_dec;

        for (uint32_t i = 0; i < input_vec.size(); i++)
        {
            // Verify decibel change is an integer.
            try
            {
                input_dec = stol(input_vec[i]);
            }
            catch (const invalid_argument &e)
            {
                throw invalid_argument("Invalid decibel argument. Expected a number. \n");
            }
            catch (const out_of_range &e)
            {
                throw out_of_range("Decibel scaler is out of range for 'std::stol'. \n");
            }

            // Verify decibel change is betwen -24 and 24.
            if (input_dec < -24 || input_dec > 24)
            {
                throw out_of_range("Decibel scaler is outside valid range. Valid rand is (-24dB, 24dB)\n");
            }

            modify.push_back(input_dec);
        }

        // Bands are set to ten frequency ranges.
        bands = get_bands();
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Flat band equalizer to apply decibel changes across entire band.
     * 
     * @param fft FFT object to apply equalizing.
     * @param fs sampling rate of signal.
     * @return vector<complex<double>> Equalized FFT values
     */
    vector<complex<double>> flat_band(fft &fft, uint32_t &fs)
    {

        vector<complex<double>> fft_data;
        vector<uint32_t> indices;

        // Convert decibel change to gain.
        double gain;

        fft_data = fft.get_data();

        // Get frequency index bound of each band.
        for (uint32_t i = 0; i < bands.size(); i++)
        {
            indices.push_back(fft.get_index(bands[i], fs));
        }

        // Apply gain across entire band.
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

    /**
     * @brief Get the preset ten frequency bands.
     * 
     * @return vector<frequency> vector of frequency bands.
     */
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

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Modifications to be made to frequency band in decibel changes.
     * 
     */
    vector<int32_t> modify;

    /**
     * @brief Preset ten frequency bands.
     * 
     */
    vector<frequency> bands;
};

//                                     End class equalize                                        //
// ============================================================================================= //

// ============================================================================================= //
//                                Begin class class periodogram                                  //

/**
 * @brief Class to compute periodogram, frequency versus power spectrum.
 * 
 */
class periodogram
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a periodogram object given an FFT object and the sampling rate of the signal.
     * 
     * @param fft_object FFT object.
     * @param f_samp sampling rate of signal.
     */
    periodogram(fft &fft_object, uint32_t &f_samp)
    {
        // Get FFT bins.
        vector<complex<double>> fft_data = fft_object.get_data();
        uint32_t freq;

        for (uint32_t i = 0; i < fft_data.size(); i++)
        {
            // Compute modulus of complex FFT bins.
            magnitude.push_back(abs(fft_data[i]));

            // Get frequency value for each FFT bins.
            freq = (uint32_t)round((double)(i * f_samp) / (double)fft_data.size());

            // For small and large values, round to audible range.
            if (freq < 20)
            {
                freq = 20;
            }
            if (freq > 20000)
            {
                freq = 20000;
            }

            frequency f_value(freq);
            freq_values.push_back(f_value);
        }
    }

    /**
     * @brief Construct a periodogram object given averaged FFT bins from windowing. 
     * 
     * @param avg_values average FFT values.
     * @param f_samp sampling rate of signal.
     */
    periodogram(vector<double> &avg_values, uint32_t &f_samp)
    {
        magnitude = avg_values;
        uint32_t freq;

        for (uint32_t i = 0; i < avg_values.size(); i++)
        {
            // Get frequency value for each FFT bins.
            freq = (uint32_t)round(((double)i * (double)f_samp) / (double)avg_values.size());
            if (freq < 20)
            {
                freq = 20;
            }
            if (freq > 20000)
            {
                freq = 20000;
            }
            frequency f_value(freq);
            freq_values.push_back(f_value);
        }
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Write power-frequency values to .csv file.
     * 
     * @param file name of output file.
     */
    void write(string &file)
    {
        ofstream output_periodogram(file);

        if (!output_periodogram.is_open())
        {
            throw ios_base::failure("Error opening output file! \n");
        };

        // Write header.
        string header = "frequency,power\n";
        output_periodogram << header;

        // Write power-frequency values.
        for (uint32_t i = 0; i < magnitude.size(); i++)
        {
            output_periodogram << freq_values[i].get_frequency() << "," << magnitude[i] << "\n";
        }
    }

    /**
     * @brief Get the magnitude object that contains the power values.
     * 
     * @return vector<double> 
     */
    vector<double> get_magnitude()
    {
        return magnitude;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Power spectrum values.
     * 
     */
    vector<double> magnitude;

    /**
     * @brief Frequency values.
     * 
     */
    vector<frequency> freq_values;
};

//                                     End class periodogram                                     //
// ============================================================================================= //

// ============================================================================================= //
//                                   Begin class process_audio                                   //
/**
 * @brief Class to perform all audio processing effects. 
 * 
 */
class process_audio
{
public:
    // ===========
    // Constructor
    // ===========

    /**
     * @brief Construct a new process audio object for filtering and equalizing.
     * 
     * @param input wav_file object.
     * @param type type of audio process to perform.
     * @param f1 first frequency object to correspond to audio process type (low, high, band).
     * @param f2 second frequency object to correspond to audio process type (band).
     * @param e equalize object to correspond to equalize type.
     */
    process_audio(wav_file &input, string &type, frequency &f1, frequency &f2, equalize &e)
    {
        process(input, type, f1, f2, e);
    }

    /**
     * @brief Construct a new process audio object for mixing effects.
     * 
     * @param input_one first wav_file object to apply effect.
     * @param input_two second wav_file object to apply effect.
     * @param type mixing effect type (add or overlap).
     */
    process_audio(wav_file &input_one, wav_file &input_two, string &type)
    {
        // Check for the same number of channels.
        if (input_one.get_channel() != input_two.get_channel())
        {
            throw invalid_argument("Input files have different number of channels. Cannot mix files.\n");
        }

        // Check for the same signal sampling rate.
        if (input_one.get_samp() != input_two.get_samp())
        {
            cout << "Input files have different sampling rates. Using larger sampling rate. \n";

            if (input_two.get_samp() > input_one.get_samp())
            {
                uint32_t f = input_two.get_samp();
                input_one.update_fsamp(f);
            }
        }

        cout << "Mixing input files...\n";
        if (type == "add")
        {
            add(input_one, input_two);
        }

        if (type == "overlap")
        {
            overlap(input_one, input_two);
        }
    }

    // =======================
    // Public member functions
    // =======================

    /**
     * @brief Mixing function to concatenate second audio file to first audio file.
     * 
     * @param input_one first wav_file object.
     * @param input_two second wav_file object.
     */
    void add(wav_file &input_one, wav_file &input_two)
    {
        vector<complex<double>> data_one;
        vector<complex<double>> data_two;

        data_one = input_one.get_data();
        data_two = input_two.get_data();

        // Total data size as the sum of the two.
        uint32_t s = input_one.orig_size() + input_two.orig_size();

        // 16-bit size to not read in padding values.
        uint32_t input_one_bitsize = input_one.orig_size() / 2;
        uint32_t input_two_bitsize = input_two.orig_size() / 2;

        // Update WAV file header with total data size.
        input_one.update_size(s);

        output_data = data_one;

        // Concatenate both data vectors.
        output_data.insert(output_data.begin() + input_one_bitsize, data_two.begin(), data_two.begin() + input_two_bitsize);
    }

    /**
     * @brief Mixing function to overlap two audio signals.
     * 
     * @param input_one first wav_file object.
     * @param input_two second wav_file object.
     */
    void overlap(wav_file &input_one, wav_file &input_two)
    {
        vector<complex<double>> data_one;
        vector<complex<double>> data_two;

        uint32_t m;                    // Data size.
        complex<double> max_value = 0; // Max amplitude in new signal.

        data_one = input_one.get_data();
        data_two = input_two.get_data();

        // Compute 16-bit data size to remove padding.
        uint32_t input_one_bitsize = input_one.orig_size() / 2;
        uint32_t input_two_bitsize = input_two.orig_size() / 2;

        // If first file has longer signal, then sum two signal and keep the remaining longer signal.
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

                // Keep track of maximum amplitude.
                if (real(output_data[i]) > real(max_value))
                {
                    max_value = output_data[i];
                }
            }
        }
        else
        {
            // Update first wav_file header with larger signal size.
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

        // If maximum amplitude is larger than 2^16/2, thant scale all amplitude values.
        if (real(max_value) > 32768)
        {
            double scale = real(max_value) / 32768;
            for (uint32_t i = 0; i < output_data.size(); i++)
            {
                output_data[i] = round(real(output_data[i]) / scale) + 0i;
            }
        }
    }

    /**
     * @brief Processing function to perform FFT, filtering/equalizing and invese FFT.
     * 
     * @param input wav_file object.
     * @param type type of audio process to perform.
     * @param f1 first frequency object to correspond to audio process type (low, high, band).
     * @param f2 second frequency object to correspond to audio process type (band).
     * @param e equalize object to correspond to equalize type.
     */
    void process(wav_file &input, string &type, frequency &f1, frequency &f2, equalize &e)
    {
        cout << "Processing WAV file...\n";

        vector<complex<double>> input_data;       // Input signal data.
        vector<complex<double>> fft_data;         // Data after FFT is applied.
        vector<complex<double>> ifft_data;        // Data after IFFT is applied.
        vector<complex<double>> transformed_data; // Data after audio processing.
        uint32_t sampling_rate;                   // Sampling rate of signal.

        // Get input data, size, channel and sampling rate.
        input_data = input.get_data();
        uint64_t input_size = input_data.size();
        uint16_t channel = input.get_channel();
        sampling_rate = input.get_samp();

        if (channel == 1)
        {
            string periodogram_name = "periodogram.csv";

            // If input size is small, no windowing is required.
            if (input_size <= 262144)
            {
                // Perform FFT on entire signal.
                fft fft_all(input_data);
                fft_data = fft_all.get_data();

                // Compute periodogram and write output.
                periodogram power(fft_all, sampling_rate);
                power.write(periodogram_name);

                // Process signal by process type.
                transformed_data = process_type(fft_all, type, f1, f2, e, sampling_rate);

                // Compute IFFT on transformed signal.
                i_fft ifft_all(transformed_data);

                // Output results.
                output_data = ifft_all.get_data();
            }
            else
            {
                // If data size is large, apply windowing to prevent running out of memory.
                output_data = window(input_size, input_data, type, f1, f2, e, sampling_rate, periodogram_name);
            }
        }
        else // If channel is stereo, apply process of left and right channels separately.
        {
            cout << "Processing for left and right channels...\n";

            vector<complex<double>> left;              // Left channel.
            vector<complex<double>> right;             // Right channel.
            vector<complex<double>> fft_left_data;     // Left channel FFT data.
            vector<complex<double>> fft_right_data;    // Right channel FFT data.
            vector<complex<double>> transformed_left;  // Processed left channel data.
            vector<complex<double>> transformed_right; // Processed right channel data.
            vector<complex<double>> output_left;       // Left channel output.
            vector<complex<double>> output_right;      // Right channel output.
            vector<double> freq_spec_left;             // Frequency-power data for left channel.
            vector<double> freq_spec_right;            // Frequency-power data for right channel.
            vector<double> average;                    // Average frequency-power data.
            uint64_t channel_size = input_size / 2;    // Size of each channel.

            // Split data into both channels.
            for (uint32_t j = 0; j < channel_size; j++)
            {
                left.push_back(input_data[j * 2]);
                right.push_back(input_data[(j * 2) + 1]);
            }
            // If input size is small, no windowing is required.
            if (channel_size <= 262144)
            {
                // Perform FFT on each channel.
                fft fft_left(left);
                fft fft_right(right);
                fft_left_data = fft_left.get_data();
                fft_right_data = fft_right.get_data();

                // Compute periodogram on each channel and write output.
                periodogram power_left(fft_left, sampling_rate);
                periodogram power_right(fft_right, sampling_rate);
                string periodogram_name = "periodogram_leftchannel.csv";
                power_left.write(periodogram_name);
                periodogram_name = "periodogram_rightchannel.csv";
                power_right.write(periodogram_name);

                // Process signal by process type.
                transformed_left = process_type(fft_left, type, f1, f2, e, sampling_rate);
                transformed_right = process_type(fft_right, type, f1, f2, e, sampling_rate);

                // Compute IFFT on transformed signal.
                i_fft ifft_left(transformed_left);
                i_fft ifft_right(transformed_right);

                // Output results.
                output_left = ifft_left.get_data();
                output_right = ifft_right.get_data();
            }
            else
            {

                // If data size is large, apply windowing to prevent running out of memory.
                string periodogram_name = "periodogram_leftchannel.csv";
                cout << "Left channel\n";
                output_left = window(channel_size, left, type, f1, f2, e, sampling_rate, periodogram_name);

                periodogram_name = "periodogram_rightchannel.csv";
                cout << "Right channel\n";
                output_right = window(channel_size, right, type, f1, f2, e, sampling_rate, periodogram_name);
            }

            // Reconstruct signal from left and right channels.
            for (uint32_t j = 0; j < channel_size; j++)
            {
                output_data.push_back(output_left[j]);
                output_data.push_back(output_right[j]);
            }
        }
    }

    /**
     * @brief Function to perform processing when signal is too large. Processing is performed by divinding up the signal into windows.
     * 
     * @param input_size signal size.
     * @param input_data signal.
     * @param type type of audio process to perform.
     * @param f1 first frequency object to correspond to audio process type (low, high, band).
     * @param f2 second frequency object to correspond to audio process type (band).
     * @param e equalize object to correspond to equalize type.
     * @param sampling_rate sampling rate of signal.
     * @param periodogram_name name of output periodogram file.
     * @return vector<complex<double>> processed data.
     */
    vector<complex<double>> window(uint64_t &input_size, vector<complex<double>> &input_data, string &type, frequency &f1, frequency &f2, equalize &e, uint32_t &sampling_rate, string &periodogram_name)
    {
        vector<double> ramp;                      // Linear ramp from 0 to 1 with length of window size
        vector<complex<double>> fft_data;         // Data after FFT is performed.
        vector<complex<double>> ifft_data;        // Data after IFFT is performed.
        vector<complex<double>> transformed_data; // Processed data.
        vector<complex<double>> final;            // Final process data.
        vector<double> freq_spec;                 // Frequency-Power data per window.
        vector<double> freq_spec_total;           // Total Frequency-Power data.

        uint64_t window_size = 262144;                   // Size to partion signal into.
        uint64_t num_windows = input_size / window_size; // Number of windows to process.
        uint32_t progress_counter = 1;                   // Progress counter.

        freq_spec_total.insert(freq_spec_total.begin(), window_size, 0);

        // Create ramp vectpr.
        ramp.push_back(0);

        for (uint32_t j = 1; j < (window_size - 1); j++)
        {
            ramp.push_back(((double)j) / ((double)window_size - 1));
        }
        ramp.push_back(1);

        // For each window
        for (uint32_t i = 0; i < (num_windows - 1); i++)
        {

            cout << "[";
            // Number of loops completed.
            for (uint32_t progress = 0; progress < progress_counter; progress++)
            {
                cout << "-";
            }
            // Number of loops remaining.
            for (uint32_t progress = progress_counter; progress < (num_windows - 1); progress++)
            {
                cout << " ";
            }
            // Percentage of loop completion.
            cout << "]" << round((double)progress_counter / (double)(num_windows - 1) * (double)100) << "% \r";

            // Create two chunks from input signal.
            vector<complex<double>> chunk_one(input_data.begin() + (i * window_size), input_data.begin() + ((i * window_size) + window_size));
            vector<complex<double>> chunk_two(input_data.begin() + ((i * window_size) + window_size), input_data.begin() + ((i * window_size) + (2 * window_size)));

            // Multiply chunks by positive and negative linear ramp
            for (uint32_t k = 0; k < chunk_one.size(); k++)
            {

                chunk_one[k] = chunk_one[k] * ramp[k];
                chunk_two[k] = chunk_two[k] - (chunk_two[k] * ramp[k]);
            }

            // Concatenate chunks into one vector.
            chunk_one.insert(chunk_one.end(), chunk_two.begin(), chunk_two.end());

            // FFT the entire chunk.
            fft fft_chunk_one(chunk_one);
            fft_data = fft_chunk_one.get_data();

            // Compute periodogram on FFT chunk.
            periodogram power(fft_chunk_one, sampling_rate);
            freq_spec = power.get_magnitude();

            // Sum power data from periodogram on each chunk.
            for (uint32_t j = 0; j < freq_spec_total.size(); j++)
            {
                freq_spec_total[j] = freq_spec_total[j] + freq_spec[j];
            }

            // Process signal according to type.
            transformed_data = process_type(fft_chunk_one, type, f1, f2, e, sampling_rate);

            // IFFT chunk.
            i_fft ifft_chunk_one(transformed_data);
            ifft_data = ifft_chunk_one.get_data();

            if (i == 0)
            {
                final.insert(final.begin(), input_data.size(), 0);
            }

            // Reconstruct signal from chunks.
            for (uint32_t k = 0; k < ifft_data.size(); k++)
            {

                final[(i * window_size) + k] = final[(i * window_size) + k] + ifft_data[k];
            }

            progress_counter++;
        }
        cout << "\n";

        // Divide total power spectrum by window size to get average.
        for (uint32_t j = 0; j < freq_spec_total.size(); j++)
        {
            freq_spec_total[j] = freq_spec_total[j] / (double)window_size;
        }

        // Write periodogram.
        periodogram power(freq_spec_total, sampling_rate);
        power.write(periodogram_name);
        return final;
    }

    /**
     * @brief Function to handle type of processing.
     * 
     * @param _fft FFT object.
     * @param type type of audio process to perform.
     * @param f1 first frequency object to correspond to audio process type (low, high, band).
     * @param f2 second frequency object to correspond to audio process type (band).
     * @param e equalize object to correspond to equalize type.
     * @param f_samp sampling rate of signal.
     * @return vector<complex<double>> Processed data.
     */
    vector<complex<double>> process_type(fft &_fft, string &type, frequency &f1, frequency &f2, equalize &e, uint32_t f_samp)
    {
        uint32_t index;                           // Initialize index object.
        uint32_t index_two;                       // Initialize second index object.
        vector<complex<double>> fft_data;         // Input FFT data.
        vector<complex<double>> transformed_data; // Output processed data.
        fft_data = _fft.get_data();

        if (type == "low")
        {
            // Get index from frequency.
            index = _fft.get_index(f1, f_samp);

            // Create filter.
            filter low_fil(index, fft_data);

            // Get low-pass filter output data.
            transformed_data = low_fil.low_pass();
        }
        if (type == "high")
        {
            // Get index from frequency.
            index = _fft.get_index(f1, f_samp);

            // Create filter.
            filter high_fil(index, fft_data);

            // Get high-pass filter output data.
            transformed_data = high_fil.high_pass();
        }
        if (type == "band")
        {
            // Get both indices from frequencies.
            index = _fft.get_index(f1, f_samp);
            index_two = _fft.get_index(f2, f_samp);

            // Create filter.
            filter band_fil(index, index_two, fft_data);

            // Get band-pass filter output data.
            transformed_data = band_fil.band_pass();
        }
        if (type == "equalize")
        {
            // Get equalized output data.
            transformed_data = e.flat_band(_fft, f_samp);
        }
        return transformed_data;
    }

    /**
     * @brief Get the processed FFT data.
     * 
     * @return vector<complex<double>> data after processing.
     */
    vector<complex<double>> get_data()
    {
        return output_data;
    }

private:
    // ============
    // Private data
    // ============

    /**
     * @brief Processed FFT data.
     * 
     */
    vector<complex<double>> output_data;
};

//                                     End class process_audio                                   //
// ============================================================================================= //

int main(int argc, char *argv[])
{
    // Check for valid number of arguments.
    if (argc == 5 || argc == 6 || argc == 14)
    {

        try
        {
            vector<complex<double>> input_data;  // Input signal data.
            vector<complex<double>> output_data; // Output signal data.
            string arg_two;                      // String for argument two.
            string arg_three;                    // String for argument three.

            arg_two = argv[2];

            // Check that argument two is valid .wav file name.
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
                    throw invalid_argument("Provided 4 arguments to program. Third argument is not \"low\" or \"high\" or \"add\" or \"overlap\". \n");
                }
                else if (arg_three == "low" || arg_three == "high")
                {
                    // For low-pass and high-pass filtering

                    // Validate frequency input.
                    frequency input_freq(argv[4]);

                    // Dummy equalize object.
                    equalize buffer;

                    // Process signal.
                    process_audio processor(input, arg_three, input_freq, input_freq, buffer);

                    // Get output signal.
                    output_data = processor.get_data();
                }
                else
                {
                    // For concatenating or overlapping signals

                    // Read and validate second WAV file.
                    wav_file input_two(argv[4]);

                    // Process signal.
                    process_audio mixer(input, input_two, arg_three);

                    // Get output signal.
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
                    // For band-pass filter.

                    // Validate frequency inputs.
                    frequency lower_freq(argv[4]);
                    frequency upper_freq(argv[5]);
                    uint32_t lower = lower_freq.get_frequency();
                    uint32_t upper = upper_freq.get_frequency();

                    if (upper <= lower)
                    {
                        throw invalid_argument("Second frequency argument is smaller than first frequency argument. \n");
                    }

                    // Dummy equalize object.
                    equalize buffer;

                    // Process signal.
                    process_audio processor(input, arg_three, lower_freq, upper_freq, buffer);

                    // Get output signal.
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
                    // For equalizing.

                    // Dummy frequency objects.
                    uint32_t buffer = 20;
                    frequency buffer_1(buffer);
                    frequency buffer_2(buffer);

                    // Initialize equalize object.
                    equalize modify(argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12], argv[13]);

                    // Process signal.
                    process_audio processor(input, arg_three, buffer_1, buffer_2, modify);

                    // Get output signal.
                    output_data = processor.get_data();
                }

                break;
            }

            cout << "Complete. \n";

            // Write output WAV file.
            wav_file output(input, output_data);
            output.write_wav(arg_two);
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
        cout << "3. The type of audio processing to perform. Valid arguments are: \"low\", \"high\", \"band\", \"equalize\", \"add\", \"overlap\"\n";
        cout << "4. The arguments to be passed to the type of audio process. \n";
        cout << "\t low: One frequency value in (20,20000) Hz to apply the low-pass filter. Example: 1000 \n";
        cout << "\t high: One frequency value in (20,20000) Hz to apply the high-pass filter. Example: 1000 \n";
        cout << "\t band: Two frequency values in (20,20000) Hz to apply the band-pass filter. Example: 1000 6000\n";
        cout << "\t equalize: Ten decibel values between (-24,24) to scale frequency bands. Example: 1 12 0 0 -3 8 -8 9 1 0  \n";
        cout << "\t add: The name of a WAV audio file to add to the first input file. Example: \"infile_two.wav\" \n";
        cout << "\t overlap: The name of a WAV audio file to be overlapped with the first input file. Example: \"infile_two.wav\" \n";
    }
}