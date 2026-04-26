# Empirical Q^k(chi_P) vs Hardy-Littlewood prediction

S_2 = 2.300938    (truncated, primes <= 100)
S_3 = 54.116088    (truncated, primes <= 47)

W-tricked predictions:
  W=6: S^(W)_2 = 1.0226, S^(W)_3 = 1.3362
  W=30: S^(W)_2 = 1.0069, S^(W)_3 = 1.0946
  W=210: S^(W)_2 = 1.0023, S^(W)_3 = 1.0292

=== N = 1024 ===
chi_P density p = 0.167969, expected pi(N)/N = 0.167969

  Q^2(chi_P)        = 2.1031
  Q^2(Bernoulli)    = 1.0492 ± 0.0066 (n=10)
  HL prediction S_2 = 2.3009
  ratio Q^2(chi)/S_2 = 0.9140    (1.0 = exact match)
  Q^3(chi_P)        = 35.6148
  Q^3(Bernoulli)    = 5.0946 ± 1.1735
  HL prediction S_3 = 54.1161
  ratio Q^3(chi)/S_3 = 0.6581

  Liouville:
  Q^2(L)            = 1.0017  (predicted ~1; Liouville centered is Gowers-uniform)
  Q^2(Bernoulli@1/2)= 1.0017 ± 0.0001

  W-trick W=6:
  Q^2(chi_W,b)      = 1.0169
  Q^2(Bernoulli)    = 1.0051 ± 0.0009
  HL prediction S^(W)_2 = 1.0226    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 0.9944
  Q^3(chi_W,b)      = 1.2869
  HL prediction S^(W)_3 = 1.3362
  ratio Q^3(chi_W)/S^(W)_3 = 0.9631

  W-trick W=30:
  Q^2(chi_W,b)      = 1.0068
  Q^2(Bernoulli)    = 1.0044 ± 0.0003
  HL prediction S^(W)_2 = 1.0069    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 0.9999
  Q^3(chi_W,b)      = 1.1375
  HL prediction S^(W)_3 = 1.0946
  ratio Q^3(chi_W)/S^(W)_3 = 1.0392

  W-trick W=210:
  Q^2(chi_W,b)      = 1.0054
  Q^2(Bernoulli)    = 1.0051 ± 0.0013
  HL prediction S^(W)_2 = 1.0023    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 1.0031
  Q^3(chi_W,b)      = 1.1370
  HL prediction S^(W)_3 = 1.0292
  ratio Q^3(chi_W)/S^(W)_3 = 1.1047

=== N = 4096 ===
chi_P density p = 0.137695, expected pi(N)/N = 0.137695

  Q^2(chi_P)        = 2.1316
  Q^2(Bernoulli)    = 1.0196 ± 0.0015 (n=10)
  HL prediction S_2 = 2.3009
  ratio Q^2(chi)/S_2 = 0.9264    (1.0 = exact match)
  Q^3(chi_P)        = 35.4396
  Q^3(Bernoulli)    = 3.1704 ± 0.2819
  HL prediction S_3 = 54.1161
  ratio Q^3(chi)/S_3 = 0.6549

  Liouville:
  Q^2(L)            = 1.0004  (predicted ~1; Liouville centered is Gowers-uniform)
  Q^2(Bernoulli@1/2)= 1.0004 ± 0.0000

  W-trick W=6:
  Q^2(chi_W,b)      = 1.0139
  Q^2(Bernoulli)    = 1.0021 ± 0.0001
  HL prediction S^(W)_2 = 1.0226    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 0.9915
  Q^3(chi_W,b)      = 1.2239
  HL prediction S^(W)_3 = 1.3362
  ratio Q^3(chi_W)/S^(W)_3 = 0.9159

  W-trick W=30:
  Q^2(chi_W,b)      = 1.0051
  Q^2(Bernoulli)    = 1.0016 ± 0.0001
  HL prediction S^(W)_2 = 1.0069    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 0.9982
  Q^3(chi_W,b)      = 1.0926
  HL prediction S^(W)_3 = 1.0946
  ratio Q^3(chi_W)/S^(W)_3 = 0.9981

  W-trick W=210:
  Q^2(chi_W,b)      = 1.0029
  Q^2(Bernoulli)    = 1.0018 ± 0.0002
  HL prediction S^(W)_2 = 1.0023    (-> 1 as W -> infty)
  ratio Q^2(chi_W)/S^(W)_2 = 1.0006
  Q^3(chi_W,b)      = 1.0659
  HL prediction S^(W)_3 = 1.0292
  ratio Q^3(chi_W)/S^(W)_3 = 1.0357

=== N = 16384 ===
chi_P density p = 0.115967, expected pi(N)/N = 0.115967

  Q^2(chi_P)        = 2.1460
  Q^2(Bernoulli)    = 1.0071 ± 0.0002 (n=10)
  HL prediction S_2 = 2.3009
  ratio Q^2(chi)/S_2 = 0.9327    (1.0 = exact match)

=== N = 65536 ===
chi_P density p = 0.099823, expected pi(N)/N = 0.099823

  Q^2(chi_P)        = 2.1489
  Q^2(Bernoulli)    = 1.0025 ± 0.0001 (n=10)
  HL prediction S_2 = 2.3009
  ratio Q^2(chi)/S_2 = 0.9339    (1.0 = exact match)