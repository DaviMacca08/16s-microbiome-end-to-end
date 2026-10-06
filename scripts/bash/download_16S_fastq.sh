#!/bin/bash

# ============================================================
# Download diretto fastq.gz (via ENA) per i 42 campioni 16S
# del dataset principale PRJNA759642 (Zuo et al. 2022)
# ============================================================
set -euo pipefail

OUTDIR="../fastq_16S_main_cohort"
mkdir -p "$OUTDIR"
cd "$OUTDIR"
echo "Directory di lavoro: $(pwd)"

RUNS=(
SRR15702528 SRR15702529 SRR15702530 SRR15702531 SRR15702532 SRR15702534
SRR15702535 SRR15702536 SRR15702537 SRR15702538 SRR15702540 SRR15702542
SRR15702544 SRR15702546 SRR15702548 SRR15702552 SRR15702554 SRR15702556
SRR15702558 SRR15702559 SRR15702560 SRR15702561 SRR15702562 SRR15702563
SRR15702564 SRR15702566 SRR15702567 SRR15702568 SRR15702628 SRR15702629
SRR15702630 SRR15702631 SRR15702632 SRR15702633 SRR15702634 SRR15702635
SRR15702636 SRR15702638 SRR15702639 SRR15702640 SRR15702641 SRR15702642
)

echo "Totale run da scaricare: ${#RUNS[@]}"

# 1) Interroga l'API ENA per l'intero progetto (query bulk per singolo run non affidabile)

echo "Interrogo ENA API (intero progetto PRJNA759642)..."
HTTP_CODE=$(curl -s -w "%{http_code}" -f \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA759642&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv" \
  -o run_report_full.tsv) || {
    echo "ERRORE: la chiamata a ENA è fallita (HTTP $HTTP_CODE)."
    echo "Controlla la connessione o riprova tra qualche minuto."
    exit 1
  }
echo "HTTP status: $HTTP_CODE"

if [ ! -s run_report_full.tsv ]; then
  echo "ERRORE: run_report_full.tsv è vuoto. Interrompo."
  exit 1
fi

echo "Righe totali nel progetto (114 run attesi + header): $(wc -l < run_report_full.tsv)"

# Filtra localmente solo i 42 run del dataset principale (awk, portabile su macOS/Linux)
printf '%s\n' "${RUNS[@]}" > wanted_runs.txt
head -n 1 run_report_full.tsv > run_report.tsv
awk -F'\t' 'NR==FNR{want[$1]=1; next} want[$1]' wanted_runs.txt run_report_full.tsv >> run_report.tsv

N_FOUND=$(($(wc -l < run_report.tsv) - 1))
if [ "$N_FOUND" -lt "${#RUNS[@]}" ]; then
  echo "ATTENZIONE: trovati $N_FOUND run su ${#RUNS[@]} attesi. Controlla wanted_runs.txt vs run_report_full.tsv."
fi

echo "--- Anteprima del report filtrato (prime 3 righe) ---"
head -n 3 run_report.tsv
echo "---------------------------------------------"
echo "Report filtrato: $(wc -l < run_report.tsv) righe (attese: 43, incluso header)"

# 2) Estrai i link FTP (ogni run ha 2 file, R1 e R2, separati da ';')
tail -n +2 run_report.tsv | while IFS=$'\t' read -r run ftp md5; do
  IFS=';' read -ra LINKS <<< "$ftp"
  for link in "${LINKS[@]}"; do
    echo "ftp://${link}"
  done
done > download_links.txt

echo "Link totali da scaricare: $(wc -l < download_links.txt)"

if [ ! -s download_links.txt ]; then
  echo "ERRORE: nessun link estratto da run_report.tsv. Controlla il contenuto del report qui sopra."
  exit 1
fi

# 3) Download parallelo (4 file alla volta: buon compromesso per M3/rete casalinga,
#    non sovraccarica la RAM perché ogni chiamata curl è leggera)
#    -O : salva con il nome file originale
#    -C -: riprende un download interrotto invece di ripartire da zero
#    -f : fallisce esplicitamente su errori HTTP (niente file vuoti/corrotti silenziosi)
cat download_links.txt | xargs -n 1 -P 4 -I{} curl -s -f -O -C - {}

echo "Download completato. File in: $(pwd)"
ls -lh *.fastq.gz | wc -l
