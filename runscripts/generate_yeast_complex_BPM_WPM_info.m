function generate_yeast_complex_BPM_WPM_info(phenos)

% phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};

getenv('BRIDGEPATH') 
projectdir=sprintf('project_yeast_%s_complex_t25_b50_mhygeSSI',phenos); 
generate_BPM_WPM_info(projectdir,'snp_pathway_min5_max300.mat','snpgenemapping_500bp.mat') 

cd(getenv('BRIDGEPATH')) 
projectdir=sprintf('project_yeast_%s_complex_t50_b25_mhygeSSI',phenos); 
generate_BPM_WPM_info(projectdir,'snp_pathway_min5_max300.mat','snpgenemapping_500bp.mat') 
