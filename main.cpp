#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream> //string stream
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
            current = input.get();
            wav_identifier.insert(i, current);
        }

        if (wav_identifier != "WAVE")
        {
            throw invalid_argument("Not a .WAV file.");
        }
        input.seekg(0);
        getline(input, s);
        //string header;
        //string data;
        header = s.substr(0, 44);
        cout << header << "\n";
        data = s.substr(44, string::npos);
        cout << data << "\n";

        // //header and data
        // for (uint32_t j = 0; j < 44; j++)
        // {
        //     header.push_back(input.get());
        //     cout << header[j] << "\n";
        // }
        // cout << *(int32_t *)&header[40] << "\n";

        // //might have to use getline(input,s)
        // //use substr to split string into header and data
        // //header[40] is an int8t
        // for (uint32_t j = 0; j < header[0]; j++)
        // {
        //     data.push_back(input.get());
        //     cout << data[j] << "\n";
        // }

        //close input file
        input.close();
    }

    string header;
    string data;
};

int main()
{
    try
    {
        wav_file input("piano.wav");
        cout << input.header << "\n";
    }
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
    }

    // //open input file
    // ifstream input("piano.wav", ios::binary);

    // if (!input.is_open())
    // {
    //     cout << "Error opening input file!";
    //     return -1;
    // }
    // //close input file
    // input.close();
}