source('src/data_analysis/extract_data.R')

cols_quanti <- c("Q2", "AGE10")

cols_quali <- c("TYPE", "AGRI", "CFA", "CAPBE", "ANTER", "SUPER", "INFG3",
                "IMP39", "BP3", "ERA", "MODTH", "REGETAB", "FRAETR09",
                "Q1", "Q16", 
                "Q25", "Q31", "SIXIEMEREG", "SIXIEMECATAEU", "SIXIEMESTATUTUU",
                "Q33", "Q34", "Q35NEW", "BACREG", "BACTAEU", "BACSTATUTUU",
                "Q39", "OS1", "OS3", 
                "PHD", "ETR1", "FP1",
                "Q50_13", "PER1", "CA7", "CA8", "CA13", "CA24", "")

# SIXIEMETAU -Tranche d'aire urbaine -> à voir si quali ou quanti