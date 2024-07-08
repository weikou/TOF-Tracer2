module massLibrary

import ..MasslistFunctions


FullPrimaryionslist_NH4soft =   [
        MasslistFunctions.massFromComposition(H=2, O=1)
        MasslistFunctions.massFromComposition(H=4, O=2)
        MasslistFunctions.massFromComposition(H=6, O=3)
        MasslistFunctions.massFromComposition(H=8, O=4)
        MasslistFunctions.massFromComposition(H=3, N=1)
        MasslistFunctions.massFromComposition(H=5, N=1, O=1)
        MasslistFunctions.massFromComposition(H=7, N=1, O=2)
        MasslistFunctions.massFromComposition(H=9, N=1, O=3)
        MasslistFunctions.massFromComposition(H=6, N=2)
        MasslistFunctions.massFromComposition(H=8, N=2, O=1)
        MasslistFunctions.massFromComposition(H=9, N=3)
        MasslistFunctions.massFromComposition(H=11, N=3, O=1)
    ]

CO2 = MasslistFunctions.createCompound(C=1, O=2)
NH3 = MasslistFunctions.createCompound(N=1, H=3)
PYRIDINE = MasslistFunctions.createCompound(N=1,H=5,C=5)
ACETONE = MasslistFunctions.createCompound(C=3,H=6,O=1)
ACETONE_nh4 = MasslistFunctions.createCompound(C=3,H=9,O=1,N=1)
MVK_nh4 = MasslistFunctions.createCompound(C=4,H=9,O=1,N=1)
APINENE = MasslistFunctions.createCompound(C=10, H=16)
APINENE_nh4 = MasslistFunctions.createCompound(C=10, H=19,N=1)
BCARY = MasslistFunctions.createCompound(C=15, H=24)
ISOPRENE = MasslistFunctions.createCompound(C=5, H=8)
HEXANONE =  MasslistFunctions.createCompound(C=6, H=12, O=1)
HEXANONE_nh4 =  MasslistFunctions.createCompound(C=6, H=15, O=1,N=1)
PINONALDEHYDE = MasslistFunctions.createCompound(C=10, H=16, O=2)
PINONALDEHYDEPAN = MasslistFunctions.createCompound(C=10, H=15, O=6, N=1)
PINONICACID = MasslistFunctions.createCompound(C=10, H=16, O=3) #C10H16O3
PINICACID = MasslistFunctions.createCompound(C=9, H=14, O=4) #C9H14O4
ACETIC = MasslistFunctions.createCompound(C=2, H=4, O=2)
ACETICFRAG = MasslistFunctions.createCompound(C=2, H=2, O=1)
NH3 = MasslistFunctions.createCompound(N=1, H=3)
DMA = MasslistFunctions.createCompound(C=2, H=7, N=1)
TMA = MasslistFunctions.createCompound(C=3, H=9, N=1)
TMB = MasslistFunctions.createCompound(C=9, H=12)
D5Siloxane = (371.10183)
D5Siloxane_NH4 = (388.12837)

ACETONITRILE = MasslistFunctions.createCompound(C=2, H=3, N=1)

OrgNitratNO = MasslistFunctions.createCompound(C=10, H=15, O=5, N=1)
NORPINONALDEHYDE = MasslistFunctions.createCompound(C=9, H=14, O=2)
NORPINONALDEHYDEPAN = MasslistFunctions.createCompound(C=9, H=13, O=6, N=1)
H3O = MasslistFunctions.createCompound(H=2, O=1)
H3OH2O = MasslistFunctions.createCompound(H=4, O=2)
H3OH2OH2O = MasslistFunctions.createCompound(H=6, O=3)

########## Sulphur ###############
DMS = MasslistFunctions.createCompound(C=2, H=6, S=1)
DMSO = MasslistFunctions.createCompound(C=2, H=6, S=1, O=1)
DMSO2 = MasslistFunctions.createCompound(C=2, H=6, S=1, O=2)
MSIA = MasslistFunctions.createCompound(C=1, H=4, S=1, O=2)
MSA = MasslistFunctions.createCompound(C=1, H=4, S=1, O=3)

########## BCARY Specific ########
#Fast
C15H22O2 = MasslistFunctions.createCompound(C=15, H=22, O=2)
C15H24O2 = MasslistFunctions.createCompound(C=15, H=24, O=2)
C15H24O3 = MasslistFunctions.createCompound(C=15, H=24, O=3)
C15H26O4 = MasslistFunctions.createCompound(C=15, H=26, O=4)
#Slow / Sticky
C14H20O3 = MasslistFunctions.createCompound(C=14, H=20, O=3)
C13H20O4 = MasslistFunctions.createCompound(C=13, H=20, O=4)
C14H22O4 = MasslistFunctions.createCompound(C=14, H=22, O=4)
C15H22O4 = MasslistFunctions.createCompound(C=15, H=22, O=4)
C15H24O4 = MasslistFunctions.createCompound(C=15, H=24, O=4)
C14H25NO4 = MasslistFunctions.createCompound(C=14, H=25, O=4, N=1)
C15H29NO6 = MasslistFunctions.createCompound(C=15, H=29, O=6, N=1)
#Nitrates
C15H23NO4 = MasslistFunctions.createCompound(C=15, H=23, O=4, N=1)
C13H19NO6 = MasslistFunctions.createCompound(C=13, H=19, O=6, N=1)
C15H25NO4 = MasslistFunctions.createCompound(C=15, H=25, O=4, N=1)
C15H21NO5 = MasslistFunctions.createCompound(C=15, H=21, O=5, N=1)
C15H23NO5 = MasslistFunctions.createCompound(C=15, H=23, O=5, N=1)
C15H23NO6 = MasslistFunctions.createCompound(C=15, H=23, O=6, N=1)
C15H17NO7 = MasslistFunctions.createCompound(C=15, H=17, O=7, N=1)

###### END BCARY ###############


###### NAPHTHA #################
NAPHTHA = MasslistFunctions.createCompound(C=10,H=8)

###### gas standards ###########

# CLOUD STD2 masses
CLOUD_STD2_masses = Dict(
	"Acetone" => ([MasslistFunctions.massFromComposition(C=3, H=6, O=1),
				MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1),
				MasslistFunctions.massFromComposition(C=3, H=8, O=2),
				MasslistFunctions.massFromComposition(C=3, H=11, O=2, N=1)
				],1000),
	"MVK" => ([MasslistFunctions.massFromComposition(C=4, H=6, O=1),
			  	MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1),
			  	MasslistFunctions.massFromComposition(C=4, H=8, O=2),
			  	MasslistFunctions.massFromComposition(C=4, H=11, O=2, N=1)
			  	],1000),
	"Hexanone" => ([MasslistFunctions.massFromComposition(C=6, H=12, O=1),
					MasslistFunctions.massFromComposition(C=6, H=15, O=1, N=1),
					MasslistFunctions.massFromComposition(C=6, H=14, O=2),
					MasslistFunctions.massFromComposition(C=6, H=17, O=2, N=1)
					],1000),
	"Heptanone" => ([MasslistFunctions.massFromComposition(C=7, H=14, O=1),
					MasslistFunctions.massFromComposition(C=7, H=17, O=1, N=1),
					MasslistFunctions.massFromComposition(C=7, H=16, O=2),
					MasslistFunctions.massFromComposition(C=7, H=19, O=2,N=1)
					],1000),
	"Apinene" => ([MasslistFunctions.massFromComposition(C=10, H=16),
				MasslistFunctions.massFromComposition(C=10, H=19, N=1),
				MasslistFunctions.massFromComposition(C=10, H=18, O=1),
				MasslistFunctions.massFromComposition(C=10, H=21, N=1, O=1)
				],1000),
	"betaCaryophyllene" => ([MasslistFunctions.massFromComposition(C=15, H=24),
				MasslistFunctions.massFromComposition(C=15, H=27, N=1),
				MasslistFunctions.massFromComposition(C=15, H=26, O=1),
				MasslistFunctions.massFromComposition(C=15, H=29, N=1, O=1)
						],1000)
)


# green STD - STD1 (2022/2023)
CLOUD_greenSTD_masses = Dict(
	"Acetonitrile" => ([MasslistFunctions.massFromComposition(C=2, H=3, N=1),
				MasslistFunctions.massFromComposition(C=2, H=6, N=2),
				MasslistFunctions.massFromComposition(C=2, H=9, N=3),
				MasslistFunctions.massFromComposition(C=2, H=5, N=1, O=1)
				],1000),
	"Acetaldehyde" => ([MasslistFunctions.massFromComposition(C=2, H=4, O=1),
				MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=1),
				MasslistFunctions.massFromComposition(C=2, H=10, N=2, O=1),
				MasslistFunctions.massFromComposition(C=2, H=9, N=1, O=2)
				],1000),
	"Ethanol" => ([MasslistFunctions.massFromComposition(C=2, H=6, O=1),
				MasslistFunctions.massFromComposition(C=2, H=9, N=1, O=1),
				MasslistFunctions.massFromComposition(C=2, H=12, N=2, O=1),
				MasslistFunctions.massFromComposition(C=2, H=11, N=1, O=2)
				],1000),
	"Acetone" => ([MasslistFunctions.massFromComposition(C=3, H=6, O=1),
				MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1),
				MasslistFunctions.massFromComposition(C=3, H=8, O=2),
				MasslistFunctions.massFromComposition(C=3, H=12, O=1, N=2),
				MasslistFunctions.massFromComposition(C=3, H=11, O=2, N=1)
				],1000),
	"Isoprene" => ([MasslistFunctions.massFromComposition(C=5, H=8),
				MasslistFunctions.massFromComposition(C=5, H=11, N=1),
				MasslistFunctions.massFromComposition(C=5, H=13, N=1, O=1),
				MasslistFunctions.massFromComposition(C=5, H=14, N=2)
				],1000),
	"MACR" => ([MasslistFunctions.massFromComposition(C=4, H=6, O=1),
			  	MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1),
			  	MasslistFunctions.massFromComposition(C=4, H=8, O=2),
			  	MasslistFunctions.massFromComposition(C=4, H=12, O=1, N=2),
			  	MasslistFunctions.massFromComposition(C=4, H=11, O=2, N=1)
			  	],1000),
	"MEK" => ([MasslistFunctions.massFromComposition(C=4, H=8, O=1),
			  	MasslistFunctions.massFromComposition(C=4, H=11, O=1, N=1),
			  	MasslistFunctions.massFromComposition(C=4, H=10, O=2),
			  	MasslistFunctions.massFromComposition(C=4, H=14, O=1, N=2),
			  	MasslistFunctions.massFromComposition(C=4, H=13, O=2, N=1)
			  	],1000),
	"Hexanone" => ([MasslistFunctions.massFromComposition(C=6, H=12, O=1),
				MasslistFunctions.massFromComposition(C=6, H=15, O=1, N=1),
				MasslistFunctions.massFromComposition(C=6, H=14, O=2),
				MasslistFunctions.massFromComposition(C=6, H=18, O=1, N=2),
				MasslistFunctions.massFromComposition(C=6, H=17, O=2, N=1)
				],1000),
	"Apinene" => ([MasslistFunctions.massFromComposition(C=10, H=16),
				MasslistFunctions.massFromComposition(C=10, H=19, N=1),
				MasslistFunctions.massFromComposition(C=10, H=18, O=1),
				MasslistFunctions.massFromComposition(C=10, H=22, N=2),
				MasslistFunctions.massFromComposition(C=10, H=21, N=1, O=1)
				],1000),
	"TMB" => ([MasslistFunctions.massFromComposition(C=9, H=12),
				MasslistFunctions.massFromComposition(C=9, H=15, N=1),
				MasslistFunctions.massFromComposition(C=9, H=18, N=2),
				MasslistFunctions.massFromComposition(C=9, H=17, N=1, O=1)
				],1000)
)

# brown STD - STD2 (2022/2023)
CLOUD_brownSTD_masses = Dict(
	"Acetonitrile" => ([MasslistFunctions.massFromComposition(C=2, H=3, N=1),
				MasslistFunctions.massFromComposition(C=2, H=6, N=2),
				MasslistFunctions.massFromComposition(C=2, H=9, N=3),
				MasslistFunctions.massFromComposition(C=2, H=5, N=1, O=1)
				],1000),
	"Acetaldehyde" => ([MasslistFunctions.massFromComposition(C=2, H=4, O=1),
				MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=1),
				MasslistFunctions.massFromComposition(C=2, H=10, N=2, O=1),
				MasslistFunctions.massFromComposition(C=2, H=9, N=1, O=2)
				],1000),
	"Acetic Acid" => ([MasslistFunctions.massFromComposition(C=2, H=4, O=2),
				MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=2),
				MasslistFunctions.massFromComposition(C=2, H=6, O=3),
				MasslistFunctions.massFromComposition(C=2, H=10, N=2, O=2),
				MasslistFunctions.massFromComposition(C=2, H=9, N=1, O=3)
				],1000),
	"Acetone" => ([MasslistFunctions.massFromComposition(C=3, H=6, O=1),
				MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1),
				MasslistFunctions.massFromComposition(C=3, H=8, O=2),
				MasslistFunctions.massFromComposition(C=3, H=12, O=1, N=2),
				MasslistFunctions.massFromComposition(C=3, H=11, O=2, N=1)
				],1000),
	"DMS" => ([MasslistFunctions.massFromComposition(C=2, H=6, S=1),
				MasslistFunctions.massFromComposition(C=2, H=9, N=1, S=1),
				MasslistFunctions.massFromComposition(C=2, H=11, N=1, S=1, O=1),
				MasslistFunctions.massFromComposition(C=2, H=12, N=2, S=1)
				],1000),
	"MVK" => ([MasslistFunctions.massFromComposition(C=4, H=6, O=1),
			  	MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1),
			  	MasslistFunctions.massFromComposition(C=4, H=8, O=2),
			  	MasslistFunctions.massFromComposition(C=4, H=12, O=1, N=2),
			  	MasslistFunctions.massFromComposition(C=4, H=11, O=2, N=1)
			  	],1000),
	"Hexenal" => ([MasslistFunctions.massFromComposition(C=6, H=10, O=1),
				MasslistFunctions.massFromComposition(C=6, H=13, O=1, N=1),
				MasslistFunctions.massFromComposition(C=6, H=12, O=2),
				MasslistFunctions.massFromComposition(C=6, H=16, O=1, N=2),
				MasslistFunctions.massFromComposition(C=6, H=15, O=2, N=1)
				],1000),
	"Hexanone" => ([MasslistFunctions.massFromComposition(C=6, H=12, O=1),
				MasslistFunctions.massFromComposition(C=6, H=15, O=1, N=1),
				MasslistFunctions.massFromComposition(C=6, H=14, O=2),
				MasslistFunctions.massFromComposition(C=6, H=17, O=2, N=1)
				],1000),
	"Octanone" => ([MasslistFunctions.massFromComposition(C=8, H=16, O=1),
				MasslistFunctions.massFromComposition(C=8, H=19, O=1, N=1),
				MasslistFunctions.massFromComposition(C=8, H=18, O=2),
				MasslistFunctions.massFromComposition(C=8, H=22, O=1, N=2),
				MasslistFunctions.massFromComposition(C=8, H=21, O=2, N=1)
				],1000),
	"Apinene" => ([MasslistFunctions.massFromComposition(C=10, H=16),
				MasslistFunctions.massFromComposition(C=10, H=19, N=1),
				MasslistFunctions.massFromComposition(C=10, H=18, O=1),
				MasslistFunctions.massFromComposition(C=10, H=22, N=2),
				MasslistFunctions.massFromComposition(C=10, H=21, N=1, O=1)
				],750),
	"Benzene" => ([MasslistFunctions.massFromComposition(C=6, H=6),
				MasslistFunctions.massFromComposition(C=6, H=9, N=1),
				MasslistFunctions.massFromComposition(C=6, H=8, O=1),
				MasslistFunctions.massFromComposition(C=6, H=12, N=2),
				MasslistFunctions.massFromComposition(C=6, H=11, N=1, O=1)
			],1000),
	"Toluene" => ([MasslistFunctions.massFromComposition(C=7, H=8),
				MasslistFunctions.massFromComposition(C=7, H=11, N=1),
				MasslistFunctions.massFromComposition(C=7, H=10, O=1),
				MasslistFunctions.massFromComposition(C=7, H=14, N=2),
				MasslistFunctions.massFromComposition(C=7, H=13, N=1, O=1)
			],1000),
	"Xylene" => ([MasslistFunctions.massFromComposition(C=8, H=10),
				MasslistFunctions.massFromComposition(C=8, H=13, N=1),
				MasslistFunctions.massFromComposition(C=8, H=12, O=1),
				MasslistFunctions.massFromComposition(C=8, H=16, N=2),
				MasslistFunctions.massFromComposition(C=8, H=15, N=1, O=1)
			],1000)
)

end
