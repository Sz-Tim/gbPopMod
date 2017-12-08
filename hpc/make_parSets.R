# set parameter sequences
library(gbPopMod)

p <- names(set_g_p())[10:22]
p.seq <- list(K=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              pr.s=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              pr.f=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              fec=expand_cnpy(Op=c(100, 200), Cl=c(10, 40), length_out=2),
              age.f=expand_cnpy(Op=c(3, 4), Cl=c(6, 7), length_out=2),
              pr.sb=seq(0.1, 0.8, length.out=3),
              pr.est=expand_cnpy(Op=c(0.01,0.1), Cl=c(0.01, 0.1), length_out=2),
              sdd.max=seq(10, 20, length.out=3),
              sdd.rate=seq(0.05, 0.8, length.out=3),
              n.ldd=seq(1, 5, length.out=3),
              pr.eat=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              bird.hab=expand_cnpy(Op=c(0.1, 0.5), Cl=c(0.1, 0.5), length_out=2),
              pr.s.bird=seq(0.3, 0.8, length.out=3))

saveRDS(p, "hpc/p.rds")
saveRDS(p.seq, "hpc/p_seq.rds")
