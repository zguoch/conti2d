
/**
 * @file MultiProgressBar.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Definition of progress bar.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef MULTIPROGRESSBAR
#define MULTIPROGRESSBAR
#include "stdio.h"
#include <iostream>
#include <vector>
#include <sys/ioctl.h>
#include <stdio.h>
using namespace std;
#define COLOR_BAR_PURPLE 0
#define COLOR_BAR_BLUE 1
#define COLOR_BAR_GREEN 2
#define COLOR_BAR_YELLOW 3
#define COLOR_BAR_RED 4
class MultiProgressBar
{
    public:
        MultiProgressBar(vector<double> xmin,vector<double> xmax,vector<string> title);
        MultiProgressBar(double total,int color=0);
        ~MultiProgressBar();
        void Update(double current_pos=-1);
        void Update(vector<double>current_pos);
    protected:
        vector<string> m_bar_str;
        int m_length_bar;
        char m_bar_char_left;
        char m_bar_char_right;
        vector<double> m_percent;
        vector<string>m_title;
        vector<double> m_total;
        vector<double> m_current_index;
        int m_bar_number;
        vector<double> m_left,m_right;
        int m_maxLength_title;
    private:
        double m_factor;
        vector<string>m_colors;
        void init_colors();
        int m_defaultcolor;
};

#endif //MULTIPROGRESSBAR

