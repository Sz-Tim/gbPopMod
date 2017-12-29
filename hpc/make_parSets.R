# set parameter sequences
library(gbPopMod)

p <- names(set_g_p())[10:22]
p.seq <- list(K=expand_LCs(g.p=g.p, p="K", LC="all", len_out=6,
                           lc.min=c(500, 10, 80, 200, 80, 80),
                           lc.max=c(800, 20, 120, 300, 120, 120)),
              pr.s=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=3),
              pr.f=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=3),
              fec=expand_cnpy(Op=c(100, 200), Cl=c(10, 40), length_out=3),
              age.f=expand_cnpy(Op=c(3, 4), Cl=c(6, 7), length_out=3),
              pr.sb=seq(0.1, 0.8, length.out=6),
              pr.est=expand_cnpy(Op=c(0.01,0.1), Cl=c(0.01, 0.1), length_out=3),
              sdd.max=round(seq(10, 20, length.out=6)),
              sdd.rate=seq(0.05, 0.8, length.out=6),
              n.ldd=round(seq(1, 6, length.out=6)),
              pr.eat=expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=3),
              bird.hab=expand_cnpy(Op=c(0.1, 0.5), Cl=c(0.1, 0.5), length_out=3),
              pr.s.bird=seq(0.3, 0.8, length.out=6))

saveRDS(p, "hpc/p.rds")
saveRDS(p.seq, "hpc/p_seq.rds")
