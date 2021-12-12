/**
 * @file main.cpp
 * @author Jennifer Freeman (freemjc@mcmaster.ca)
 * @brief 
 * @version 0.1
 * @date 2021-12-11
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream> //string stream
#include <cmath>   //is this the right math file
#include <complex>
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

        string s;
        string wav_identifier;
        string current;
        uint64_t data_size = 0;
        uint64_t padded_data_size = 0;
        uint64_t power = 1; //what happens if power is too large
        uint64_t num_16bits = data_size / 2;

        ifstream input(input_string, ios::binary);

        if (!input.is_open())
        {
            cout << "Error opening input file!";
        };

        input.seekg(8);

        for (uint32_t i = 0; i < 4; i++)
        {
            current = (char)input.get();
            wav_identifier.insert(i, current);
        }

        if (wav_identifier != "WAVE")
        {
            throw invalid_argument("Not a .WAV file.");
        }

        input.seekg(0);
        getline(input, s);

        //header is 44 bytes
        header = s.substr(0, 44);
        cout << header << "\n";

        data_size = *(int32_t *)&header[40]; //size in bytes
        cout << data_size;

        //check if mod2, else there is a padding byte
        if (data_size % 2 != 0)
        {
            data_size = data_size + 1;
        }

        while (power < num_16bits) //what to do if num_16bits is too large
        {
            power *= 2;
        }
        cout << data_size << "\n";
        cout << num_16bits << "\n";
        cout << power << "\n";

        int16_t temp;
        cout << "char " << sizeof(char) << "\n";
        cout << "int16_t " << sizeof(int16_t) << "\n";
        input.seekg(44);
        for (uint64_t i = 0; i < num_16bits; i++) //think this is reading in things right
        {

            input.read((char *)&temp, sizeof(int16_t)); //possibly working
            data.push_back((int16_t)temp);
        }

        //difference between power and audio size
        padded_data_size = power - num_16bits;

        cout << "padding size" << padded_data_size << "\n";
        data.resize(power);

        input.close();
    }

    /**
     * @brief Size of data chunk in WAV file.
     * 
     * @return int32_t 
     */
    int32_t datasize()
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
        return new_vec; //is this bad, making a copy?
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

    //might want to add, sampling rate, data size etc.
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

        const double pi = acos(-1.0L);

        //get data from wav_file
        uint64_t N = dat.size();
        if (N == 1)
        {
            //do nothing
        }
        else
        {

            double N_2 = N / 2;
            // cout << "N_2" << N_2 << "\n";
            vector<complex<double>> even;
            vector<complex<double>> odd;

            //do I need to reserve?
            //even.reserve(N_2);
            //odd.reserve(N_2);

            for (uint64_t i = 0; i < N_2; i++)
            {
                even.push_back(dat[2 * i]);
                odd.push_back(dat[2 * i + 1]);
            }

            //recursive
            cooley_tukey(even);
            cooley_tukey(odd);

            complex<double> q = 0;
            for (uint64_t k = 0; k < N_2; k++)

            {

                q = polar(1.0, (-2.0 * pi * (double)k) / (N)) * odd[k];
                dat[k] = even[k] + q;
                dat[k + N_2] = even[k] - q;
            }

            cout << "end of fft \n";
        }

        return dat;
    }

    //might not need this function
    void write_file()
    {
        FILE *file_ptr = NULL;

        file_ptr = fopen("fft.txt", "w"); //open file for writing
        if (file_ptr == NULL)
        {
            perror("Cannot open writing file");
        }

        for (uint64_t jj = 0; jj < x.size(); jj++) //write each row of chain to file
        {
            if (fprintf(file_ptr, "%f %f\n", real(x[jj]), imag(x[jj])) < 0)
            {
                perror("Error occurred");
                fclose(file_ptr);
            }
        }

        fclose(file_ptr);
    }
    //might not need this function
    void print()
    {
        for (uint64_t i = 0; i < x.size(); i++)
        {
            double current = abs(x.at(i));
            cout << current << "\n";
        }
    }
    //might not need this function
    void min_max()
    {
        double min_val = 2000000000;
        double max_val = 0;
        for (uint64_t i = 0; i < x.size(); i++)
        {
            double current = abs(x[i]);
            if (current < min_val)
            {
                min_val = current;
            }
            else if (current > max_val)
            {
                max_val = current;
            }

            cout << current << "\n";
        }
        cout << "max val " << max_val << "\n";
        cout << "min val " << min_val << "\n";
    }

    /**
     * @brief Get the data object.
     * 
     * @return vector<complex<double>> 
     */
    vector<complex<double>> get_data()
    {
        return x; //is this bad, making a copy?
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
     * @param dft values to perform inverse FFT on.
     * @return vector<complex<double>> inverse FFT values.
     */
    vector<complex<double>> inverse(vector<complex<double>> &dft)
    {
        //complex conjugate
        for (uint64_t i = 0; i < dft.size(); i++)
        {
            dft[i] = conj(dft[i]);
        }

        //fft
        fft out(dft);

        dft = out.get_data();
        int32_t N = dft.size();

        //complex conjugate
        //divide by N
        for (uint64_t i = 0; i < dft.size(); i++)
        {
            dft[i] = conj(dft[i]) / (double)N;
        }

        return dft;
    }

    //not sure if I need this member function
    void write_file()
    {
        FILE *file_ptr = NULL;

        file_ptr = fopen("ifft.txt", "w"); //open file for writing
        if (file_ptr == NULL)
        {
            perror("Cannot open writing file");
        }

        for (uint64_t jj = 0; jj < y.size(); jj++) //write each row of chain to file
        {
            if (fprintf(file_ptr, "%f %f\n", real(y[jj]), imag(y[jj])) < 0)
            {
                perror("Error occurred");
                fclose(file_ptr);
            }
        }

        fclose(file_ptr);
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

int main()
{
    try
    {
        vector<complex<double>> input_data;
        vector<complex<double>> output_data;

        // Read and validate WAV file.
        wav_file input("440_sine.wav");

        cout << input.datasize();
        input_data = input.get_data();

        // Compute fast fourier tranform on input WAV file.
        fft test(input_data);
        cout << "end of algo \n";

        //test.min_max();
        test.write_file();

        // Compute inverse fast fourier transform on input transformed values.
        output_data = test.get_data();
        i_fft output(output_data);

        output.write_file();
    }
    //catch others, error opening file...
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
    }

    //To Do
    //1. Check exceptions fo all functions being used
    //2. Use index instead of .at() where makes sense
    //3. Restrict input file to be of certain size, or split it
    //4.
}