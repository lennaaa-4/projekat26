 Prvi dan: 6.9.2026.

 Paketi potrebni za rad su ucitani u program. Set podataka takodje je ucitan u program i zatim je matrica sa transkriptima definisana kao SeuratObject (anotacija kao metapodaci). Uradjena je kontrola kvaliteta prvo na osnovu procenta mitohondrijalne DNK i onda pomocu napravljenog ViolinPlota u kome su gledani podaci o broju transkripata, broju gena. Filtriranje nije uradjeno jer je iz VlnPlota zakljuceno da su podaci spremni za dalju obradu.

 Drugi dan 7.9.2026.

 Odradjena je logaritamska normalizacija podataka kako bi se stabilizovala varijansa u sirovoj ekspresiji. Odredjeni su najvarijabilniji geni (napravljena je i lista 10 najvarijabilnijih) i na osnovu toga uradjen je VariableFeaturesPlot na kome su oblezeni prvih 10 najvarijabilnijih gena i moja ciljena RNK. VariableFeatures kasnije su skalirani kako bi bile smanjene razliike u ekspresiji gena (tj. dovedene na uporediv nivo). Takodje je odredjena i PCA, analiza kljucnih komponenti, kojom se redukuje dimenzionalnost i odredjene su kljucne bioloske razlike (varijacije) prema kojima su geni podeljeni u PC-eve. Na osnovu PCA napravljen je DimHeatMap za 500 celija u PC_1 komponenti koji sluzi da bi se proverilo kakav je zapravo bioloski signal i da li ima ostataka nekevrste tehnickog suma. Napravljen je ElbowPlot pomocu kog je odredjeno koji ce se PC-evi koristiti u daljoj analizi. Takodje je provereno prisustvo LINC00958 u PC-evima. Nakon toga uradjeno je matematicko klasterovanje celija i uradjeno je pet DimPlotova sa razlicitim rezolucijama koji sluze za podesavanje broja klastera (grupa) koji ce biti stvoreni.

 Treci dan 8.9.2026.

 Nakon analize DimPlotova odluceno je da se pokusa da se podaci vizualizuju na druge nacine jer su prethodno napravljeni plotovi bili prilicno nepregledni. To je uradjeno na tri razlicita nacina: dva puta 3D prikazom, PCA bubble plotom sto podrazumeva vizualizaciju PC_1, PC_2 (na osama) i PC_3 (kroz velicinu tacaka) i 2D plotom na koji je nadodata i treca, z-osa za PC_3. Sprovedena je UMAP analiza i napravljen je UMAP DimPlot radi vizualizacije klastera. Napisan je izvestaj za prva tri dana.

 Cetvrti dan 9.9.2026.

 Uradjen je pseudobulk koji je podrazumevao vizualizaciju, pravljenje i definisanje pocetne matrice zatim pravljenje podgrupa koje su takodje morale biti sredjene pre upotrebe, definisanje i uredjivanje metapodatka; i nakon toga uradjena je DESeq2 analiza, medjutim nakon filtracije i zavrsetka analize celokupnu obradu prosla je samo jedna RNK (LINC00486), te je doneta odluka da se narednih dana to analizira.

 Peti dan 10.9.2026.

 Citanje radova i istrazivanje o LINC00486 koja je jedina prosla analizu diferencijalne ekspresije.

 Sesti dan 11.9.2026. 

 Ponavljanje analize diferencijalne ekspresije. Odluceno je da se dodaju i kontrolne grupe (zdrava bela masa) u analizu te uporeda izgleda kao: chronic active vs control, periplaque vs control a ne periplakna vs hronicna aktivna. Nakon dodavanja kontrola, analiza je uspesna sa 11 diferencijalno eksprimiranih gena izmedju periplakne bele mase i kontrola, dok izmedju hronicno aktivnih lezija i kontrola ima 56 diferencijalno eksprimiranih gena. Za ove analize uradjeni su VolcanoPlotovi radi vizualizacije.

 Sedmi dan 12.9.2026.

 Pisanje izvestaja i komentara za figure. Uradjena je provera dijagnostickog potencijala kao biomarkera (paket pROC) za dugu nekodirajucu RNK U91319.1 koja je bila medju diferencijalno eksprimiranim i kod periplakne i kod hronicno aktivne bele mase. Rezultati dobijeni vizualizacijom su u skladu sa ocekivanjima, tj. prate logiku rezultata iz prethodne analize diferencijalne ekspresije. Odluceno je da se ne koriste ovi rezultati i analiza u projektu jer je velika mogucnost greske zbog malog broja uzoraka.

 Osmi dan 13.9.2026.

 Odradjena je nizvodna analiza signalnih puteva i njihove aktivnosti koristeci se paketom PROGENy i nappravljeni su Pheatmapovi za vizualizaciju rezultata u formi uzorci periplakna i uzorci kontrola + svi signalni putevi i uzorci hronicno aktivna i uzorci kontrola + svi signalni putevi. Nakon toga je zbog prethodne hipoteze projekta uradjena BoxPlot vizualizacija za JAK-STAT signalni put.

 Deveti dan 14.9.2026.

Uradjena je procena dijagnosticke vrednosti kao biomarkera za signalni put EGFR, koji je potvrdjeno jedan od kljucnih "regulatora" bolesti, koristeci pakete decoupleR i pROC. Uradjen je Multiclass ROC sa sve tri grupe tkiva tj. sva tri moguca para tkiva. Rezultati su prikazani u vidu ROC kriva. Daljom obradom (proverom AUC i CI skorova) utvrdjeno je da rezultati nisu statisticki znacajni.

Deseti dan 15.9.2026.

Uradjena je pROC analiza za sve signalne puteve i prikazani su rezultati koji su podrazumevali ime signalnog puta, ime para tkiva koja su uporedjena, AUC vrednost (vrednost povrsine ispod ROC krive), donja i gornja granica CI vrednosti (stepen pouzdanosti AUC procene) kao i razlika izmedju granica. U ovim rezultatima najbolji skor imao je WNT signalni put kod uporede periplakne bele mase i kontrola (AUC: 0.917, raspon CI: 0.314). Odluceno je uraditi jos nizvodnih analiza targetirajuci specificno WNT signalni put.

Jedanaesti dan 16.9.2026.

Uzeti su rezultati analize diferencijalne ekspresije i WNT vrednosti dobijene u prethodnoj analizi. WNT rezultati filtrirani su tako da se uzimaju u obzir samo geni koji imaju pozitivan ili negativan doprinos signalnom putu. Nakon toga, rezultati spojeni su sa rezultatima analize diferencijalne ekspresije i uradjen je grafik. Ovo je ponovljeno za sva tri para. Grafici oznacavaju znacaj gena tj. koji su geni najodgovorniji za promene WNT skora. Potom je istrazen znacaj WNT signalnog puta u multiploj sklerozi i trazena je ideja za preformulisanje poente rada.
