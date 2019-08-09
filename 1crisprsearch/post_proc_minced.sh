# minced.out.gff files require some post processing, because
# voglio eliminare tutte le lines che sono ###
# il post  processing di eliminare l'estensione .txt dai ,minced.out.txt l'ho gia fatto
#e è un comando, non lo scrivo
cd /shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out
for file in *.minced.out.gff; do  echo $file; sed '/^#/ d' <$file> ${file}_cut; mv ${file}_cut $file; done

# -d xke se no ti chiede ogni volta di sostituire il file
