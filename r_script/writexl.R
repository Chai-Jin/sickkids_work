install.packages("writexl")
library(writexl)

Idents(G4_whole_final) <- factor(Idents(sc_obj), levels = c("add sample here"))
Idents(G4_whole_final) <- factor(Idents(G4_whole), levels = c("G4_A_P_0527", "G4_A_R_0527", "G4_D_P_1303", "G4_D_R_1303", "G4_E_P_2126", "G4_E_R_2126"))

Idents(G4_whole_final) <- G4_whole_final$orig.ident