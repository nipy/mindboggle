#! /bin/csh -fe
# 
# apply_colormap.csh - process the files in the folder
#

set DATA = ( \
input/images/PeriSylvian.eps \
input/images/Dorsal-Ventral_surface1.eps \
input/images/Mid-Lateral_surface3.eps \
input/images/TableOfAbbrv.eps \
)

while ( $#DATA != 0 )
        set F = $DATA[1]
	set R = $F:r
	set E = $F:e

	echo python recolor_eps.py --mapFile output/parcLabels.xml \
		${R}_prepared.$E ${R}_recolored.$E 
	python recolor_eps.py --mapFile output/parcLabels.xml \
		${R}_prepared.$E ${R}_recolored.$E 

        shift DATA
end
