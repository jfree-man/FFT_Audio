#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream> //string stream
#include <cmath>   //is this the right math file
#include <complex>
using namespace std;

class wav_file
{
public:
    //maybe input should be a string stream, instead of file
    //constructor
    wav_file(const string &input_string)
    {
        read_wav(input_string);
    }

    void read_wav(const string &input_string)
    {
        //open input file
        ifstream input(input_string, ios::binary);

        if (!input.is_open())
        {
            cout << "Error opening input file!";
        }
        string s;
        string wav_identifier;
        string current;
        //getline(input, s);
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

        //does string stream help?
        //string data_string;
        //header is 44 bytes
        header = s.substr(0, 44);
        cout << header << "\n";
        uint64_t data_size = 0;
        uint64_t padded_data_size = 0;
        uint64_t power = 1; //what happens if power is too large

        data_size = *(int32_t *)&header[40]; //size in bytes
        cout << data_size;
        //does data size have to be mod 2
        //check if mod2, else there is a padding byte
        if (data_size % 2 != 0)
        {
            data_size = data_size + 1;
        }
        uint64_t num_16bits = data_size / 2;
        while (power < num_16bits) //what to do if num_16bits is too large
        {
            power *= 2;
        }
        cout << data_size << "\n";
        cout << num_16bits << "\n";
        cout << power << "\n";

        //char temp;
        int16_t temp;
        cout << "char " << sizeof(char) << "\n";
        cout << "int16_t " << sizeof(int16_t) << "\n";
        input.seekg(44);
        for (uint64_t i = 0; i < num_16bits; i++) //think this is reading in things right
        {

            input.read((char *)&temp, sizeof(int16_t)); //possibly working
            data.push_back((int16_t)temp);
            //data.push_back((double)temp);

            //cout << "data: " << data[i] << "\n";
        }

        //difference between power and audio size
        padded_data_size = power - num_16bits;

        cout << "padding size" << padded_data_size << "\n";
        data.resize(power);

        //data.insert(data.end(), padded_data_size, 0); //data size is a power of 2, catch exceptions here?

        //close input file
        cout << data.size() << "\n";
        cout << data.at(262140) << "\n";

        input.close();
    }
    // int16_t assign(const uint32_t &i) //dumb name, should rename
    // {
    //     return data.at(i);
    // }
    int32_t datasize()
    {
        return data.size();
    }

    vector<complex<double>> get_data()
    {
        vector<complex<double>> new_vec(data.begin(), data.end());
        return new_vec; //is this bad, making a copy?
    }

private:
    string header;
    //is it unsigned?
    vector<int16_t> data;
    //vector<complex<double>> data;
    //might want to add, sampling rate, data size etc.
};

class fft
{
public:
    fft(wav_file &audio)
    {
        vector<complex<double>> dat = audio.get_data();
        x = cooley_tukey(dat);
    }

    vector<complex<double>> cooley_tukey(vector<complex<double>> &dat)
    {
        //vector<complex<double>> fft_vec;
        const double pi = acos(-1.0L);
        //get data size from wav_file
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

            //complex<double> even_dft = 0;
            //complex<double> odd_dft = 0;

            //const complex<double> imag = 1.0i;

            //even.reserve(N_2);
            //odd.reserve(N_2);
            // cout << even.max_size() << "\n";
            // cout << odd.max_size() << "\n";
            //do I need to -1, I don't think so
            // cout << "N is " << N << "\n";
            //double N = 166454;

            for (uint64_t i = 0; i < N_2; i++)
            {
                even.push_back(dat[2 * i]);
                odd.push_back(dat[2 * i + 1]);
                //cout << "i" << i << "\n";
            }

            cooley_tukey(even);
            cooley_tukey(odd);

            // cout << "N is " << N << "\n";
            //cout << "even is " << even.size() << "\n";
            //cout << "odd is " << odd.size() << "\n";

            complex<double> q = 0;
            for (uint64_t k = 0; k < N_2; k++)

            {
                //new code here
                q = polar(1.0, (-2.0 * pi * (double)k) / (N)) * odd[k];
                //cout << even[k] + q << "\n";
                dat.at(k) = even.at(k) + q;
                dat[k + N_2] = even[k] - q;

                q = 0; //probably don't need to do this
            }

            cout << "end of fft \n";
        }
        // cout << "size of x" << x.size() << "\n";

        // if (x.size() > 0)
        // {
        //     cout << "x at 0" << x.at(0) << "\n";
        // }

        return dat;
    }
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
    void print()
    {
        for (uint64_t i = 0; i < x.size(); i++)
        {
            double current = abs(x.at(i));
            cout << current << "\n";
        }
    }
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

private:
    vector<complex<double>> x;
};

int main()
{
    try
    {
        cout << "size of complex double" << sizeof(complex<double>);
        wav_file input("440_sine.wav");
        //cout << input.header << "\n"; can;t access member function

        cout << input.datasize();
        vector<complex<double>> dat = input.get_data();
        fft test(input);
        cout << "end of algo \n";

        //test.min_max();
        test.write_file();

        //frequncy resolution is sample freq/(num FFT values)
    }
    //catch others, error opening file...
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
    }

    //To Do
    //1. Check exceptions fo all functions being used
    //2. Use index instead of .at() where makes sense
    //3. change input of fft to vector
    //4. inverse fft
}