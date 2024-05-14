#ifndef UTIL_H
#define UTIL_H

#include <cstdlib>
#include <cstdio>
#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>

using namespace std;

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

inline char complement(char base) {
    switch(base) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with(string const & value,  string const & starting)
{
    if (starting.size() > value.size()) return false;
    return  equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with(string const & value,  string const & ending)
{
	if (ending.size() > value.size()) return false;
	return  equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str)
{
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos)
    {
        return string{""};
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos)
    {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

[[maybe_unused]] inline int split(const string& str, vector<string>& ret_, string sep = ",")
{
    if (str.empty())
    {
        return 0;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos)
    {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos)
        {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        }
        else
        {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline string replace(const string& str, const string& src, const string& dest)
{
    string ret;

    string::size_type pos_begin = 0;
    string::size_type pos       = str.find(src);
    while (pos != string::npos)
    {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos       = str.find(src, pos_begin);
    }
    if (pos_begin < str.length())
    {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline int find_with_right_pos(const string & str, const string & pattern, int start = 0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline size_t substr_count(const string & str, const string & substr)
{
    size_t count{0};
    int pos{-1};
    while ((pos = str.find(substr, ++pos)) != string::npos)
        count++;
    return count;
}

inline string basename(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if(pos == filename.length()-1)
        return ""; // a bad filename
    else
        return filename.substr(pos+1, filename.length() - pos - 1);
}

inline string dirname(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else
        return filename.substr(0, pos+1);
}

inline string join_path(const string& dirname, const string& basename){
    if (dirname[dirname.length()-1] == '/') {
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

//Check if a string is a file or directory
inline bool file_exists(const  string& s)
{
    bool exists = false;
    if (!s.empty()) {
        struct stat status{};
        int result = stat( s.c_str(), &status );
        if (result == 0) {
            exists = true;
        }
    }
    return exists;
}


// check if a string is a directory
inline bool is_directory(const  string& path)
{
    bool isdir = false;
    struct stat status{};
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat( path.c_str(), &status );
    if ( status.st_mode &  S_IFDIR  ) {
        isdir = true;
    }
// #endif
    return isdir;
}

inline void check_file_valid(const  string& s) {
    if(!file_exists(s)){
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if(is_directory(s)){
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string
inline  string str_keep_alpha(const  string& s)
{
    string new_str;
    for (char it : s) {
        if (isalpha(it) ) {
            new_str += it;
        }
    }
    return new_str;
}


// Remove invalid sequence characters from a string
inline void str_keep_valid_sequence(  string& s, bool forceUpperCase = false)
{
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for( size_t it =0; it < s.size(); it++) {
        char c = s[it];
        if (forceUpperCase && c>='a' && c<='z') {
            c -= case_gap;
        }
        if (isalpha(c) || c == '-' || c == '*' ) {
            s[total] = c;
            total ++;
        }
    }
    s.resize(total);
}

inline void str2upper(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline size_t hamming_distance(const string& str1, const string& str2) {
    size_t diff = 0;
    auto len1 = str1.length();
    auto len2 = str2.length();
    for (auto i = 0; i < len1 && i < len2; ++i) {
        if (str1[i] != str2[i])
            diff++;
    }
    diff += len1 > len2 ? len1 - len2 : len2 - len1;
    return diff;
}

inline tuple<bool, bool, size_t, size_t> clipMatchConvert(const string& cigar_string) {
    size_t idx_start{0}, font_consume_len{0}, back_consume_len{0}, cigar_count{0};
    for (auto ci = 0; ci < cigar_string.size(); ++ci) {
        char cigar_op = cigar_string[ci];
        if (isalpha(cigar_op)) {
            if (cigar_op == 'S') {
                long op_len = stol(cigar_string.substr(idx_start, ci - idx_start));
                if (ci == cigar_string.size() - 1) {
                    back_consume_len += op_len;
                } else {
                    font_consume_len += op_len;
                }
            } else if (cigar_op == 'M' || cigar_op == 'I') {
                long op_len = stol(cigar_string.substr(idx_start, ci - idx_start));
                font_consume_len += op_len;
                back_consume_len += op_len;
            } else if (unlikely(cigar_op != 'D')) {
                return make_tuple(false, false, 0, 0);
            }
            ++cigar_count;
            idx_start = ci + 1;
        }
    }
    if (font_consume_len > back_consume_len + 2) {
        return make_tuple(true, false, back_consume_len, cigar_count);
    } else if (font_consume_len + 2 < back_consume_len) {
        return make_tuple(false, true, font_consume_len, cigar_count);
    } else {
        return make_tuple(false, false, 0, cigar_count);
    }
}

inline int64_t genome_cover_from_cigar(const string& cigar_string) {
    int64_t genome_size{0}, idx_start{0};
    for (auto ci = 0; ci < cigar_string.size(); ++ci)
        if (isalpha(cigar_string[ci])) {
            if (cigar_string[ci] == 'M' || cigar_string[ci] == 'D')
                genome_size += stol(cigar_string.substr(idx_start, ci - idx_start));
            idx_start = ci + 1;
        }
    return genome_size;
}

inline char num2qual(int num) {
    if (num > 127 - 33)
        num = 127 - 33;
    if (num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

#endif /* UTIL_H */
