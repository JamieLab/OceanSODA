library(seacarb)

dic = c(1020.0, 1050.0, 1080.0)
at = c(1000.0, 1010.0, 1020.0)
output = carb(15, at, dic, S=35, T=20.0, Patm=1.0, Pt=0, k1k2='x')
