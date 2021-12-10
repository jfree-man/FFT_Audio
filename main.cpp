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
        int32_t data_size = 0;
        int32_t padded_data_size = 0;
        int32_t power = 1;

        data_size = *(int32_t *)&header[40]; //size in bytes
        cout << data_size;
        //does data size have to be mod 2
        //check if mod2, else there is a padding byte
        if (data_size % 2 != 0)
        {
            data_size = data_size + 1;
        }
        int32_t num_16bits = data_size / 2;
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
        for (uint32_t i = 0; i < num_16bits; i++) //think this is reading in things right
        {

            input.read((char *)&temp, sizeof(int16_t)); //possibly working
            data.push_back((int16_t)temp);

            //cout << "data: " << data[i] << "\n";
        }

        //close input file
        cout << data.size() << "\n";

        input.close();
    }
    int16_t assign(const uint32_t &i)
    {
        return data.at(i);
    }

private:
    string header;
    //is it unsigned?
    vector<int16_t> data;
    //might want to add, sampling rate, data size etc.
};

class fft
{
public:
    fft(wav_file &audio)
    {
        cooley_tukey(audio);
    }

    void cooley_tukey(wav_file &audio)
    {
        vector<complex<double>> even;
        vector<complex<double>> odd;
        const double pi = acos(-1.0L);
        //const complex<double> imag = 1.0i;

        double N_2 = 166454 / 2; //do I need to -1, I don't think so
        double N = 166454;

        complex<double> even_dft = 0;
        complex<double> odd_dft = 0;

        for (double i = 0; i < N_2; i++)
        {
            for (double j = 0; j < N_2; j++) //sum over N/2-1, to capture last value N/2-1
            {

                //polar(r,theta)
                //m=j, k=i
                even_dft = even_dft + ((double)audio.assign(2.0 * j) * polar(1.0, (-2.0 * pi * j * i) / (N_2)));
                odd_dft = odd_dft + ((double)audio.assign((2.0 * j) + 1) * polar(1.0, (-2.0 * pi * j * (i + 1)) / (N_2)));
            }
            even.push_back(even_dft);
            odd.push_back(odd_dft);
            // cout << "even dft" << even_dft << "\n";
            // cout << "odd dft" << odd_dft << "\n";
            even_dft = 0;
            odd_dft = 0;
        }

        cout << "end of fft \n";
        complex<double> q = 0;

        for (double k = 0; k < N_2; k++)
        {
            q = polar(1.0, (-2.0 * pi * k) / (N)) * odd.at(k);
            cout << even.at(k) + q << "\n";
            x.insert(x.begin() + k, even.at(k) + q);
            x.insert(x.begin() + (k + N_2), even.at(k) - q);

            cout << "x at " << k << "is" << x.at(k) << "\n";
            cout << "x at " << k + N_2 << "is" << x.at(k + N_2) << "\n";
            cout << "capacity: " << x.capacity() << "\n";
        }

        cout << "end of fft \n";
    }

private:
    vector<complex<double>> x;
};

int main()
{
    try
    {
        wav_file input("sine.wav");
        //cout << input.header << "\n"; can;t access member function
        fft test(input);
    }
    //catch others, error opening file...
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
    }

    //To Do
    //1. Check exceptions fo all functions being used
    //2. Use .at() instead of index[] for all vectors, and catch exceptions
    //3. Add padding bytes to wav_file.data
    //4. Add recursion to fft
}