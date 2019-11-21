# compile
gcc-4.9.4/bin/g++ -std=c++11  -I /usr/local/genome/boost-1.67.0/include -o sim PTV_count_simulations.cpp

# 1000 times simulations for each parameter setting
var_ary=(1.00e+01 3.16e+01 1.00e+02 3.16e+02 1.00e+03)
mut_ary=(3.16e-7 1.00e-7 3.16e-8 1.00e-8 3.16e-9 1.00e-9 3.16e-10)
for var in ${var_ary[@]}; do for mut in ${mut_ary[@]}; do
    ./sim $mut $var 0 0 1000 0 100000 > s0_h0_mut${mut}_var${var}.txt; done; done

# draw figure (fig.s1) in R
source("/betelgeuse01/analysis/hamanaka/function/ggplot_setting.R")
library(tidyverse)
library(cowplot)
var_ary=c("1.00e+01","3.16e+01","1.00e+02","3.16e+02","1.00e+03")
mut_ary=c("3.16e-7","1.00e-7","3.16e-8","1.00e-8","3.16e-9","1.00e-9","3.16e-10")
n.file=length(var_ary)*length(mut_ary)
df = wrap(mut_ary,var_ary,n.file)
df2 = df %>% nest(-mut) %>% mutate(fig=map2(data,as.list(mut),~wrap_plot(df=.x,title=.y)))
plot_grid(df2$fig[[2]],df2$fig[[3]],df2$fig[[4]],df2$fig[[5]],df2$fig[[6]],df2$fig[[7]])
ggsave("1030_grid.png",width=5.25,height=3.55)

## functions ##
wrap = function(mut_ary,var_ary,n.file){
    obs_path_base="/betelgeuse01/analysis/hamanaka/new_model/sim_tmp/s0_h0_pen1_fp0_mut"
    df = data.frame(matrix(rep(NA,n.file*3),nrow=n.file))
    colnames(df) = c("mut","var","obs.n")
    i=1
    for(mut in mut_ary){
        for(var in var_ary){
            obs_path = paste0(obs_path_base,mut,"_var",var,"_100.txt")
        	obs = read_tsv(obs_path,col_names=c("n.variant","sum.freq"))
        	obs.n.mean = sum(obs$n.variant) / nrow(obs) 
            df[i,] = c(mut,var,obs.n.mean)
        	i=i+1
        }
    }    
    mutate_at(df,vars(mut,var,obs.n.mean),funs(as.numeric(.))) %>% 
        mutate(total_rate=mut*var,corr=obs.n/total_rate) %>% 
        as_tibble %>%
        return()
}

wrap_plot = function(df,title){
    ggplot(df,aes(x=log10(total_rate),y=log10(obs.n))) +
        geom_point() +
        geom_line(size=1,alpha=0.2) +
    	theme(legend.position=c(0.8,0.35),legend.title=element_text(size=4),legend.text=element_text(size=4)) +
        labs(x="Per-gene total mutability",y="Per-gene N of variant",fill="Per-site mutability") +
        geom_abline(intercept=6.75,slope=1,linetype='dashed',alpha=0.2) + 
        xlim(-9.5,-3.5) + ylim(-2.5,3) + ggtitle(title) + g
        #if(save){ggsave(out,width=1.75,height=1.75)}
}
