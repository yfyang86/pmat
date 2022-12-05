#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <map>
#include <random>
#include "utils.h"

bool inx_tie_incr(inx_tie x, inx_tie y){
    return (x.value < y.value);
}

bool inx_tie_decr(inx_tie x, inx_tie y){
    return (y.value < x.value);
}

static std::map<std::string, int> solver_sam_load_config_setting(){
/*
a: 0.67672
b: -0.03785
c: 0.
d: 0.
*/
    std::map<std::string, int> s_config;
    s_config["a"] = 0;
    s_config["b"] = 1;
    s_config["c"] = 2;
    s_config["d"] = 3;
    s_config["useBartlett"]=4;
    s_config["useFoldedNormal"]=5;
    return s_config;
}


void solver_sam_load_config(double *  config){
  std::string line;
  std::ifstream config_file("config.txt");
  std::string option;
  double option_val;
  static std::map<std::string, int> s_congfig = solver_sam_load_config_setting();
  std::map<std::string, int>::iterator it;

  for (size_t j=0; j<4; j++) config[j] = 0.;

  if (!config_file.is_open()) {
      std::cerr << "Could not open the file - '"
           << "config.txt" << "'" << std::endl;
      exit(EXIT_FAILURE);
  }else {
      while ( getline (config_file,option, ':') && getline (config_file,line) ){
                std::istringstream(line) >> option_val;
                it = s_congfig.find(option);
                if (it!=s_congfig.end()) config[it->second] = option_val;
    }
  config_file.close();
  }
}

/** BH method Randomize break the tie
 * @brief
 * By applying the above formulation, the FDR controlling MCP introduced in 
 * Benjamini and Hochberg (1995) can be rephrased into a FDR $p$-value 
 * adjustment. Benjamini and Hochberg's MCP is as follows: compare each 
 * ordered $p$-value $p_{(i)}$ with $\alpha \cdot i / m$. 
 * Let $k=\max \left\{i: p_{(i)} \leqslant \alpha \cdot i / m\right\}$, 
 * if such $k$ exists reject $\mathrm{H}_{(1)}^0, \ldots, \mathrm{H}_{(k)}^0$.
 * A $q$ sized $\mathrm{MCP}$ rejects an hypotheses $\mathrm{H}_{(i)}^0$ 
 * if for some $k, k=i, \ldots, m, p_{(k)} \cdot m / k \leqslant \alpha$, 
 * thus we define the BH FDR $p$-value adjustment as
 * $$
 * Q_{\text {adj }}^{\mathrm{BH}}\left(p_{(i)}\right)
 *     =
 *         \min _{i \leqslant k}\left\{p_{(k)} \cdot m / k\right\}
 * $$
 * @example{
 *  std::vector<double> ex;
 *  for (auto v : padj_bh(ex)){
 *    std::cout<<v<<"\t";
 *   }
 * }
*/
std::vector<double> padj_bh(std::vector<double> & x, bool rngtie){
    order_tie tmp_tie;
    int cnt = 1;
    std::map<int, order_tie> xp;
    std::vector<inx_tie> yinx;
    std::vector<inx_tie> roinx;
    inx_tie xxxx;
    int m = ((int) x.size());
    std::vector<int> xorder;
    std::vector<int> xorderinc;
    std::vector<double> y = x;


    for (auto xui:x){
        xxxx.inx = cnt - 1;
        xxxx.value = xui;
        yinx.push_back(xxxx);
        cnt++;
    }    

    std::sort(yinx.begin(), yinx.end(), inx_tie_decr);//decreasing

    cnt = 1;
    for (auto xui:yinx){
        xxxx.inx = cnt - 1;
        xxxx.value = (double) xui.inx;
        roinx.push_back(xxxx);
        cnt++;
    }    

    std::sort(roinx.begin(), roinx.end(), inx_tie_incr);//increasing

    for (int i = 0; i< m; ++i){// n log(n)
        xorder.push_back(i);
        xorderinc.push_back(i);
    }

    for (int i = 0; i< m; ++i){
        xorder[i]  = yinx[i].inx; // n
        xorderinc[i] = roinx[i].inx;
    }


    for (int i = 0; i< m; ++i){
        x[i] = yinx[i].value * m /(m-i+.0);
    }


    for (int i = 0; i< m; ++i){
      if (i+1 < m && x[i+1]>x[i]){
        x[i+1] = x[i];
      }
    }



    for (int i = 0; i< m; ++i){
        y[i] = x[xorderinc[i]];
                y[i] = (
          (y[i] > 1. ?  1. : y[i])
        );
    }

    return y;
}
