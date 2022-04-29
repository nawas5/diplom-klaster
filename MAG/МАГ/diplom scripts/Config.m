function [ config ] = Config( )
    
    config.c = 299792458;
    config.T = 1;
    config.sigma_t = 1*1e-9;
    config.sigma_n = config.sigma_t * config.c;
    config.sigma_phi = 0.1;
    config.sigma_ksi_v = 0.01;
    config.sigma_ksi_phi = 1;
    
    config.posts = [0 0 0; 10 0 0; 0 10 0; 10 10 0]';
    config.posts = [-5076.25700228714 12312.9342219411 -8241.14728795100 -1.11097161563963;-11487.2777145567 3509.94201846127 5345.24999149516 -0.445039107080202;160.134190817221 106.173298698501 188.845909395533 124.399999890768];
    config.posts = [-5e3 12e3 -8e3 0; -11e3 3.5e3 5e3 0; 10 10 10 10];
    config.posts_number = size(config.posts, 2);
    config.hei = 1*10e3;
    config.sigma_ksi = 100;
    config.sigma_h = 0;
    config.max_coord = 100e3;
    config.max_V = 600;
    config.max_acc = 100;


end

