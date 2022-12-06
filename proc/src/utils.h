#ifndef MAIN_UTILS_H
#define MAIN_UTILS_H

#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <map>
#include <random>
#include <map>
#include <algorithm>
#include <vector>

static std::map<std::string, int> solver_sam_load_config_setting();
void solver_sam_load_config(double * conifg);

struct  order_tie
{
   double value;
   std::vector<int> loc;
   int cnt;
};

struct inx_tie{
    double value;
    int inx;
};

bool inx_tie_incr(inx_tie x, inx_tie y);
bool inx_tie_decr(inx_tie x, inx_tie y);

std::vector<double> padj_bh(std::vector<double> & x);

#endif