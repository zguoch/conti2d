/**
 * @file MultiProgressBar.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Implementation of progress bar.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "MultiProgressBar.h"
#include <iomanip>
#include "Conti2D.h"

void MultiProgressBar::init_colors()
{
    m_colors.push_back("\033[35m"); //Purple
    m_colors.push_back("\033[34m"); //blue
    m_colors.push_back("\033[32m"); //green
    m_colors.push_back("\033[33m"); //yellow
    m_colors.push_back("\033[31m"); //red
}
MultiProgressBar::~MultiProgressBar()
{
}
MultiProgressBar::MultiProgressBar(double total, int color) : m_bar_char_left('#'),
                                                              m_bar_char_right('-'),
                                                              m_defaultcolor(color)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    m_length_bar = w.ws_col - 35;
    init_colors();
    m_total.push_back(total);
    cout << endl;
    string bar_str = "";
    for (int i = 0; i < m_length_bar; i++)
        bar_str += m_bar_char_right;
    m_bar_str.push_back(bar_str);
    m_current_index.push_back(0);
    m_percent.push_back(0);
    m_left.push_back(0);
    m_right.push_back(total);
    m_title.push_back("");

    m_factor = m_length_bar / 100.0;
}
MultiProgressBar::MultiProgressBar(vector<double> left, vector<double> right, vector<string> title) : m_bar_char_left('#'),
                                                                                                      m_bar_char_right('-'),
                                                                                                      m_title(title)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    m_length_bar = w.ws_col - 35;

    init_colors();
    if (right.size() != left.size())
    {
        cout << "The length of left and right in MultiProgressBar are not the same" << endl;
        exit(0);
    }
    m_maxLength_title = 0;
    for (int i = 0; i < left.size(); i++)
    {
        double total = fabs(left[i] - right[i]);
        m_total.push_back(total);
        cout << endl;
        string bar_str = "";
        for (int i = 0; i < m_length_bar; i++)
            bar_str += m_bar_char_right;
        m_bar_str.push_back(bar_str);
        m_current_index.push_back(0);
        m_percent.push_back(0);
        m_left.push_back(left[i]);
        m_right.push_back(right[i]);
        if (m_maxLength_title < m_title[i].size())
            m_maxLength_title = m_title[i].size();
    }

    m_factor = m_length_bar / 100.0;
}

void MultiProgressBar::Update(double current_pos)
{
    for (int k = 0; k < m_total.size(); k++)
    {
        m_current_index[k] = (current_pos < 0 ? m_current_index[k] + 1 : current_pos);
    }
    Update(m_current_index);
}
void MultiProgressBar::Update(vector<double> current_pos)
{
    if (!((current_pos.size() == m_total.size()) && (m_total.size() == m_title.size())))
    {
        cout << "The size of current_pos, m_total, m_title have different size in MultiProgressBar" << endl;
        exit(0);
    }
    for (int k = 0; k < m_total.size(); k++)
    {
        cout << "\r";
        if (k == 0)
            MOVEUP(m_total.size());
        // m_current_index[k]++;
        m_percent[k] = (fabs((current_pos[k]) - m_left[k]) / m_total[k]) * 100;
        //how many are completed
        int pos = int(m_percent[k] * m_factor);
        for (int i = 0; i <= pos; i++)
            m_bar_str[k][i] = m_bar_char_left;
        for (int i = pos + 1; i < m_length_bar; i++)
            m_bar_str[k][i] = m_bar_char_right;
        if (m_title[k] == "")
        {
            cout << "[" << m_colors[(k + m_defaultcolor) % m_colors.size()]
                 << m_bar_str[k]
                 << COLOR_DEFALUT << "]"
                 << m_colors[(k + m_defaultcolor) % m_colors.size()] << setw(3) << right << int(m_percent[k]) << "%"
                 << COLOR_DEFALUT << endl;
        }
        else
        {
            cout << "[" << m_colors[(k + m_defaultcolor) % m_colors.size()]
                 << m_bar_str[k]
                 << COLOR_DEFALUT << "]"
                 << m_colors[(k + m_defaultcolor) % m_colors.size()] << setw(3) << right << int(m_percent[k])
                 << "% " << COLOR_DEFALUT
                 << "[" << setw(m_maxLength_title) << left << m_title[k]
                 << "] [" << m_colors[(k + m_defaultcolor) % m_colors.size()] << setw(10) << left << current_pos[k] << COLOR_DEFALUT
                 << "]"
                 << endl;
        }
    }
}