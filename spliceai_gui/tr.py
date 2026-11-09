"""Türkçe texts for i18n.t(), keyed by the English text in the code.

Owner's rules (2026-11-08): short and clear; "ekzon", never "ekson"; terms and
proper names stay as they are (masked, raw, build, live, precomputed,
in-frame, frameshift, WT, cDNA, Het/Hom/Hemi, SpliceAI, MANE Select, SnpEff,
HGVS, VCF, FASTA); results-table column headers and their tooltips stay
English. Glossary: splice site = splice bölgesi, sample = örnek, acceptor /
donor = akseptör / donör, cryptic = kriptik, exon skipping = ekzon atlama,
intron retention = intron tutulumu, pseudoexon = psödoekzon, WT predicted =
WT tahmini, candidate = aday, event = olay, zygosity = zigosite.
tests/test_i18n.py checks every t("...") in the code has an entry here.
"""

TR = {
    # --- align.py
    'The text holds {n} FASTA sequences; paste one.':
        'Metinde {n} FASTA dizisi var; bir tane yapıştırın.',
    'There is text before the FASTA header line (>...).':
        'FASTA başlık satırından (>...) önce metin var.',
    'The sequence contains characters other than A, C, G, T and N: {shown}.':
        'Dizi A, C, G, T ve N dışında karakterler içeriyor: {shown}.',
    'The sequence has {n} bases; at least {least} are needed to place it reliably on the reference.':
        'Dizi {n} baz; referansa güvenilir şekilde yerleştirmek için en az {least} baz gerekir.',
    "{text} isn't a region like chr7:117,480,025-117,668,665.":
        '{text}, chr7:117,480,025-117,668,665 gibi bir bölge değil.',
    '{text}: the start has to be at least 1 and not after the end.':
        '{text}: başlangıç en az 1 olmalı ve bitişten sonra olmamalı.',
    '{text}: a position is at least 1.':
        "{text}: bir konum en az 1'dir.",
    '{text} is neither a region like chr7:117,480,025-117,668,665 nor a position like chr7:117,559,590.':
        '{text}, ne chr7:117,480,025-117,668,665 gibi bir bölge ne de chr7:117,559,590 gibi bir konum.',
    "The window ends at {end:,}, past the chromosome's end ({length:,}).":
        'Pencere {end:,} konumunda bitiyor, kromozomun sonunu ({length:,}) geçiyor.',
    "The sequence doesn't align to the region.":
        'Dizi bölgeye hizalanmıyor.',
    'Identity {pct:.2f}% over {n:,} aligned columns ({method}).':
        '{n:,} hizalı sütunda benzerlik %{pct:.2f} ({method}).',
    'Covered reference: {start:,}..{end:,} ({n:,} bases; bases {q0:,}..{q1:,} of the sequence).':
        'Kapsanan referans: {start:,}..{end:,} ({n:,} baz; dizinin {q0:,}..{q1:,} bazları).',
    "The sequence doesn't cover the whole gene: {parts} come from the reference and are treated as normal.":
        'Dizi genin tamamını kapsamıyor: {parts} referanstan alınır ve normal kabul edilir.',
    '{a} and {b}':
        '{a} ve {b}',
    'The sequence fits {build} better ({pct:.2f}% identity there): check the genome build.':
        "Dizi {build} ile daha iyi uyuşuyor (orada benzerlik %{pct:.2f}): genom build'ini kontrol edin.",
    "SpliceAI's gene interval ({s0:,}..{s1:,}) and the canonical transcript ({c0:,}..{c1:,}) differ; the analysis covers both.":
        "SpliceAI'nin gen aralığı ({s0:,}..{s1:,}) ile kanonik transkript ({c0:,}..{c1:,}) farklı; analiz ikisini de kapsar.",
    '{n} difference(s) found.':
        '{n} fark bulundu.',
    '{kind_a} at {a:,} and {kind_b} at {b:,} overlap; choose one.':
        '{a:,} konumundaki {kind_a} ile {b:,} konumundaki {kind_b} örtüşüyor; birini seçin.',
    'The {kind} at {pos:,} reaches outside the analysed region ({start:,}..{end:,}).':
        '{pos:,} konumundaki {kind}, analiz edilen bölgenin ({start:,}..{end:,}) dışına taşıyor.',
    'The {kind} at {start:,} touches the {name} ({pos:,}).':
        '{start:,} konumundaki {kind}, {name} bölgesine ({pos:,}) dokunuyor.',
    # --- cdna.py
    'start shifted {n:+d}':
        'başlangıç {n:+d} kaymış',
    'end shifted {n:+d}':
        'bitiş {n:+d} kaymış',
    'identical':
        'aynı',
    '{status} (one exon with reference exon {exons})':
        '{status} ({exons} numaralı referans ekzonuyla tek ekzon)',
    'extra (novel)':
        'fazladan (yeni)',
    # --- cdna_analysis.py
    "There is no curated RefSeq transcript (NM_) at {gene}'s locus in the {build} table, so no reference cDNA can be built, and there is no transcript start, end and ATG to anchor predicted transcripts. Nothing can be predicted for this gene.":
        '{build} tablosunda {gene} lokusunda küratörlü bir RefSeq transkripti (NM_) bulunmamaktadır; bu nedenle referans cDNA oluşturulamaz ve tahmin edilen transkriptlerin sabitleneceği transkript başlangıcı, bitişi ve ATG belirlenemez. Bu gen için tahmin yapılamaz.',
    '{text} Its MANE Select transcript {mane} exists, but the UCSC {build} table lacks it here.':
        '{text} MANE Select transkripti {mane} var, ama UCSC {build} tablosunda burada yok.',
    'No MANE Select transcript is available for this gene on {build}. The reference cDNA and the anchors (transcript start, end and ATG) are taken from {accession}.':
        'Bu gen için {build} üzerinde MANE Select transkripti yok. Referans cDNA ve çapalar (transkript başlangıcı, bitişi ve ATG) {accession} transkriptinden alınır.',
    'Found within {span:,} bases of {position:,}: {where}.':
        '{position:,} konumunun {span:,} bazı içinde bulundu: {where}.',
    'Not found within {span:,} bases of {position:,}; found in the whole gene: {where}.':
        '{position:,} konumunun {span:,} bazı içinde bulunamadı; genin tamamında bulundu: {where}.',
    'Not found within {span:,} bases of {position:,}, nor in the whole gene.':
        '{position:,} konumunun {span:,} bazı içinde de, genin tamamında da bulunamadı.',
    'No event given: one per line, e.g. "deletion chr17:17228330-17228462".':
        'Olay verilmedi: her satıra bir tane, ör. "deletion chr17:17228330-17228462".',
    'Built from {n} event(s) given as coordinates, with no alignment: the sample is the reference {region} with them applied, and the reference bases were checked.':
        'Koordinat olarak verilen {n} olaydan, hizalama yapılmadan oluşturuldu: örnek, bu olayların uygulandığı {region} referansıdır ve referans bazları kontrol edildi.',
    '{n} difference(s).':
        '{n} fark.',
    'Reference values (database or computed live)':
        'Referans değerleri (veritabanından ya da live hesaplanarak)',
    'Computing the reference live':
        'Referans live hesaplanıyor',
    'Sample haplotype: predicting near the edits':
        'Örnek haplotipi: değişikliklerin yakını tahmin ediliyor',
    'Sample haplotype':
        'Örnek haplotipi',
    'Transcripts, cDNA and protein':
        'Transkriptler, cDNA ve protein',
    'Sample rank 1 leads rank 2 by {gap:.3f} in relative score ({a:.3f} vs {b:.3f}).':
        "Örneğin 1. sıradaki adayı 2. sıradakinin göreli skorda {gap:.3f} önünde ({a:.3f}'e karşı {b:.3f}).",
    'kept in one exon by the WT prediction (retained)':
        'WT tahmini tek ekzonda tutuyor (tutulum)',
    'an exon of {n:,} nt in the WT prediction only':
        "yalnızca WT tahmininde olan {n:,} nt'lik bir ekzon",
    'acceptor of exon {n}':
        'ekzon {n} akseptörü',
    'donor of exon {n}':
        'ekzon {n} donörü',
    'before the transcript':
        'transkriptten önce',
    'after the transcript':
        'transkriptten sonra',
    'Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction matches the reference for exons {exons} ({n} of {n}).':
        '{pos:,} konumundaki {kind} çevresinde ({size:,} bp, {where}): WT tahmini {exons} ekzonlarında referansla örtüşüyor ({n}/{n}).',
    'exon {k} ({status})':
        'ekzon {k} ({status})',
    'intron {k} ({status})':
        'intron {k} ({status})',
    'Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction differs from the reference for {which}, and matches it for exons {same} ({n_same} of {n} exons identical).':
        '{pos:,} konumundaki {kind} çevresinde ({size:,} bp, {where}): WT tahmini şuralarda referanstan farklı: {which}; {same} ekzonlarında örtüşüyor ({n} ekzondan {n_same} tanesi aynı).',
    'Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction differs from the reference for {which} (0 of {n} exons identical).':
        '{pos:,} konumundaki {kind} çevresinde ({size:,} bp, {where}): WT tahmini şuralarda referanstan farklı: {which} ({n} ekzondan 0 tanesi aynı).',
    'Low wild-type concordance: only {same} of {n} reference exons are predicted identically for the reference sequence. Predictions for the sample are less reliable for this gene.':
        'Düşük wild type uyumu: referans dizisi için {n} referans ekzondan yalnızca {same} tanesi doğru tahmin edilmektedir. Bu gende örnek tahminleri daha az güvenilirdir.',
    '{text} Interpret the sample prediction there with caution.':
        '{text} Oradaki örnek tahminini dikkatle yorumlayın.',
    'Weak evidence: every splice site where sample candidate {rank} differs from the WT prediction changed by less than {delta:g}.':
        "Zayıf kanıt: {rank}. örnek adayının WT tahmininden farklı olduğu her splice bölgesi {delta:g}'den daha az değişti.",
    '{utr} in exon {n}':
        'ekzon {n} içinde {utr}',
    'intron {n}':
        'intron {n}',
    'upstream of the transcript':
        'transkriptin önünde (upstream)',
    'downstream of the transcript':
        'transkriptin ardında (downstream)',
    'between {a} and {b}':
        '{a} ile {b} arasında',
    'spanning exons {a} to {b}':
        'ekzon {a}–{b} boyunca',
    'spanning {a} to {b}':
        '{a} – {b} boyunca',
    'spanning {a} to {b} (exons {n})':
        '{a} – {b} boyunca (ekzon {n})',
    'spanning {a} to {b} (exon {n})':
        '{a} – {b} boyunca (ekzon {n})',
    '{site}: at it':
        '{site}: tam üzerinde',
    '{site}: {d:,} bp':
        '{site}: {d:,} bp',
    'inverted copy of {region} (in place)':
        '{region} bölgesinin ters kopyası (yerinde)',
    'inverted copy of {region}':
        '{region} bölgesinin ters kopyası',
    'templated from {region}':
        '{region} bölgesinden kopyalanmış',
    'bases {a}-{b}: {what}':
        'baz {a}-{b}: {what}',
    'tandem copy of {region}':
        '{region} bölgesinin tandem kopyası',
    'too short to trace (under {n} bases)':
        'kaynağı izlenemeyecek kadar kısa ({n} bazdan az)',
    'novel (not found in the gene region)':
        'yeni (gen bölgesinde bulunamadı)',
    '{region}; whole transcript covered':
        '{region}; transkriptin tamamı kapsanıyor',
    '{gap} (exons {n})':
        '{gap} (ekzonlar {n})',
    '{gap} (exon {n})':
        '{gap} (ekzon {n})',
    '{region}; missing {a} and {b}':
        '{region}; eksik: {a} ve {b}',
    '{region}; missing {a}':
        '{region}; eksik: {a}',
    'the sample was built from the reference with the events given.':
        'örnek, verilen olaylarla referanstan oluşturuldu.',
    "the sequence was given in transcript orientation (5'→3'), as in Ensembl/NCBI gene exports, and is used as such.":
        "dizi, Ensembl/NCBI gen dışa aktarımlarındaki gibi transkript yönünde (5'→3') verildi ve öyle kullanılıyor.",
    'the sequence was given as the reference (+) strand, as in genome browsers, and is used as such.':
        'dizi, genom tarayıcılarındaki gibi referans (+) zinciri olarak verildi ve öyle kullanılıyor.',
    "the sequence was given as its reverse complement (3'→5') and is used turned around.":
        "dizi ters-tamamlayıcı (3'→5') olarak verildi ve çevrilerek kullanılıyor.",
    '{gene} is on the {strand} strand; {how}':
        '{gene} {strand} zincirinde; {how}',
    'Identity':
        'Benzerlik',
    'Covered region':
        'Kapsanan bölge',
    'Orientation':
        'Yön',
    'Differences':
        'Farklar',
    # --- cdna_help.py
    "In short: the tab shows, one below the other, the reference exon boundaries, the transcript SpliceAI predicts for the normal (WT) sequence and the transcripts it predicts for the sample, so that the effect of a change can be read from their differences.\n\nThis tab predicts how a sample's DNA sequence for a single gene may be spliced, and the cDNA and protein that would result. The sample sequence is entered by pasting it, by loading it from a FASTA file, or by building it from events given as coordinates; the tab then shows three kinds of transcript, one below the other.\n\nThe reference transcript is the gene's annotated canonical transcript (MANE Select on hg38). The WT predicted transcript is predicted from the reference genome sequence alone: SpliceAI assigns each base a probability of being a splice acceptor and a splice donor, and the strongest sites are chained into transcripts. Only the transcript start, the transcript end and the ATG are taken from the annotation; all exon boundaries in between are predicted. The transcripts predicted for the sample are obtained by applying the same method to the sample sequence; up to five candidates are shown, ranked from the most likely.\n\nComparing the WT prediction with the reference shows how well the method reproduces the known transcript of this gene. Where the two agree, a difference between the sample prediction and the WT prediction is more likely to reflect the sample's change than a limitation of the method. Where the WT prediction already differs from the reference, a predicted change in that region should be interpreted with more caution.":
        "Kısaca: sekme, referans ekzon sınırlarını, SpliceAI'ın normal (WT) dizi için öngördüğü transkripti ve örnek dizi için öngördüğü transkriptleri alt alta gösterir; bir değişikliğin etkisi bu üçü arasındaki farklardan okunur.\n\nBu sekme, tek bir gene ait örnek DNA dizisinin nasıl splice edilebileceğini ve bunun sonucunda oluşacak cDNA ile proteini öngörür. Örnek dizi, yapıştırılarak, bir FASTA dosyasından yüklenerek ya da koordinatla verilen olaylardan oluşturularak girilir; sekme ardından üç tür transkripti alt alta gösterir.\n\nReferans transkript, genin anotasyondaki kanonik transkriptidir (hg38'de MANE Select). WT tahmini, yalnızca referans genom dizisinden öngörülen transkripttir: SpliceAI her baz için splice akseptörü ve splice donörü olma olasılığını hesaplar; en güçlü bölgeler birleştirilerek transkriptler oluşturulur. Anotasyondan yalnızca transkriptin başlangıcı, bitişi ve ATG alınır; aradaki tüm ekzon sınırları tahmin edilir. Örnek dizi için öngörülen transkriptler, aynı yöntem örnek diziye uygulanarak elde edilir; en olasıdan başlayarak en fazla beş aday gösterilir.\n\nWT tahmininin referansla karşılaştırılması, yöntemin bu genin bilinen transkriptini ne ölçüde doğru ürettiğini gösterir. İkisinin örtüştüğü bölgelerde, örnek tahmini ile WT tahmini arasındaki bir farkın yöntemin sınırlılığından çok örnekteki değişiklikten kaynaklanması daha olasıdır. WT tahmininin referanstan zaten farklı olduğu bölgelerde ise öngörülen değişiklikler daha temkinli yorumlanmalıdır.",
    'Site threshold':
        'Splice bölgesi eşiği',
    'The lowest SpliceAI probability at which a base can be used as a splice site when transcripts are built.':
        'Transkriptler oluşturulurken bir bazın splice bölgesi olarak kabul edilebilmesi için gereken en düşük SpliceAI olasılığıdır.',
    'A higher value keeps only strong sites: fewer and clearer candidates, at the risk of missing a weak but real exon (if one of its sites falls below the threshold, the exon can only be predicted as skipped). A lower value admits weaker sites: more candidates, with more noise.':
        'Daha yüksek bir değer yalnızca güçlü bölgeleri dikkate alır; daha az ve daha belirgin aday elde edilir, ancak zayıf fakat gerçek bir ekzon gözden kaçabilir (bölgelerinden biri eşiğin altında kalırsa ekzon yalnızca atlanmış olarak öngörülebilir). Daha düşük bir değer zayıf bölgeleri de dahil eder; aday sayısı artar, ancak gürültü de artar.',
    'Minimum internal exon':
        'En kısa iç ekzon',
    'The shortest internal exon a predicted transcript may contain. Micro-exons (internal exons of 3 to about 30 nucleotides) are common in brain and muscle transcripts; exon 3 of TP53, for example, is 22 nucleotides long. An exon shorter than this value cannot be predicted, even when SpliceAI scores both of its sites strongly. The default of 3 nt keeps micro-exons within reach: compared with 25 nt, it predicts 686 more exons correctly across all hg38 genes while adding 88 incorrect ones. SpliceAI tends to score micro-exons lower than longer exons, so not all of them are found. The first and last exons are not subject to this limit.':
        "Tahmin edilen bir transkriptte bulunabilecek en kısa iç ekzondur. Mikro-ekzonlar (3 ile yaklaşık 30 nükleotit arasındaki iç ekzonlar) beyin ve kas transkriptlerinde sık görülür; örneğin TP53'ün 3. ekzonu 22 nükleotit uzunluğundadır. Bu değerden kısa bir ekzon, SpliceAI her iki bölgesini de yüksek skorlasa bile tahmin edilemez. Varsayılan 3 nt değeri mikro-ekzonların da bulunabilmesini sağlar: 25 nt ile karşılaştırıldığında tüm hg38 genlerinde 686 ekzon daha doğru tahmin edilirken 88 hatalı ekzon eklenir. SpliceAI mikro-ekzonları uzun ekzonlara göre daha düşük skorlama eğiliminde olduğundan bunların tamamı bulunamayabilir. İlk ve son ekzon bu sınırlamaya tabi değildir.",
    'A higher value removes very short exons, including real micro-exons and a few incorrect ones. Values below 3 nt are of no use, as no annotated hg38 exon is that short.':
        "Daha yüksek bir değer, gerçek mikro-ekzonlar ve az sayıda hatalı ekzon dahil olmak üzere çok kısa ekzonları dışarıda bırakır. Hiçbir anotasyonlu hg38 ekzonu bu kadar kısa olmadığından 3 nt'nin altındaki değerler bir yarar sağlamaz.",
    'Minimum intron':
        'En kısa intron',
    'The shortest intron a predicted transcript may have.':
        'Tahmin edilen bir transkriptin sahip olabileceği en kısa intron.',
    'A higher value excludes short introns and, with them, some cryptic sites close to an exon. A lower value allows very short stretches to be spliced out, which rarely occurs in reality.':
        'Daha yüksek bir değer kısa intronları ve bunlarla birlikte ekzona yakın bazı kriptik bölgeleri dışarıda bırakır. Daha düşük bir değer çok kısa bölgelerin de splice edilerek çıkarılmasına izin verir; bu durum gerçekte nadiren görülür.',
    'Apply by default from':
        'Varsayılan uygulama sınırı',
    'Differences of this size or larger found by Check compatibility, as well as all inversions and duplications, are selected under Apply. Smaller differences are not selected, because they often represent sequencing noise or known polymorphisms (these are scored in the VCF tab).':
        'Uyumluluğu kontrol et ile bulunan bu büyüklükteki veya daha büyük farklar ile tüm inversiyon ve duplikasyonlar Uygula sütununda seçili gelir. Daha küçük farklar seçilmez; bunlar çoğunlukla dizileme gürültüsü veya bilinen polimorfizmlerdir (bu tür varyantlar VCF sekmesinde skorlanır).',
    'A higher value selects fewer differences by default; a lower value selects more. Each difference can still be selected or cleared individually.':
        'Daha yüksek bir değer varsayılan olarak daha az farkı, daha düşük bir değer daha fazla farkı seçer. Her fark ayrıca tek tek seçilebilir veya seçimi kaldırılabilir.',
    'Transcripts per sequence (Advanced)':
        'Dizi başına transkript (Gelişmiş)',
    'The number of candidate transcripts shown for the WT prediction and for the sample.':
        'WT tahmini ve örnek için gösterilecek aday transkript sayısıdır.',
    'A higher value also shows less likely alternatives; the relative scores are then shared among more transcripts, so each becomes smaller. A lower value shows only the most likely few.':
        'Daha yüksek bir değer olasılığı daha düşük alternatifleri de gösterir; bu durumda göreli skorlar daha fazla transkript arasında paylaşıldığından her birinin değeri küçülür. Daha düşük bir değer yalnızca en olası birkaç transkripti gösterir.',
    'Weak evidence below delta (Advanced)':
        'Zayıf kanıt delta sınırı (Gelişmiş)',
    'A sample candidate is marked "weak evidence" when every splice site at which it differs from the WT prediction changed by less than this value. Jaganathan et al. (Cell, 2019) used delta-score cut-offs of 0.2 (high recall: most splice-altering variants pass, with more false positives), 0.5 (their recommended balance) and 0.8 (high precision: fewer false positives, more variants missed).':
        'WT tahmininden farklı olduğu tüm splice bölgelerindeki değişim bu değerin altında kalan örnek adayları "zayıf kanıt" olarak işaretlenir. Jaganathan ve ark. (Cell, 2019) delta skoru için 0,2 (yüksek duyarlılık: splicing\'i etkileyen varyantların çoğu yakalanır, yanlış pozitifler artar), 0,5 (önerilen denge) ve 0,8 (yüksek kesinlik: yanlış pozitifler azalır, gözden kaçan varyantlar artar) eşiklerini kullanmıştır.',
    'A higher value marks more candidates as weak evidence; a lower value marks fewer. The ranking is not affected.':
        'Daha yüksek bir değer daha fazla adayı, daha düşük bir değer daha az adayı zayıf kanıt olarak işaretler. Sıralama bundan etkilenmez.',
    'Random stand-ins (Advanced)':
        'Rastgele kontrol dizileri (Gelişmiş)',
    'The number of random replacements the sequence-specificity control uses for each change.':
        'Dizi özgüllüğü kontrolünde her değişiklik için kullanılan rastgele kontrol dizisi sayısıdır.',
    'More stand-ins give a more stable comparison but take longer, since each requires one more prediction (about a third of a second for a typical gene). Fewer stand-ins are faster but less stable.':
        'Daha fazla kontrol dizisi daha kararlı bir karşılaştırma sağlar, ancak her biri ek bir tahmin gerektirdiğinden süre uzar (tipik bir gen için yaklaşık üçte bir saniye). Daha az kontrol dizisi daha hızlıdır, ancak sonuç daha az kararlıdır.',
    'Relative score':
        'Göreli skor',
    "Each transcript's score is the log-likelihood of its splice sites under SpliceAI's probabilities: each candidate site the transcript uses contributes log(p), and each candidate site it ignores contributes log(1 - p). Skipping a strong exon or retaining an intron is therefore penalised by the strong sites it ignores. The relative score converts the scores of the transcripts shown into shares that sum to 1; it is a ranking, not an estimate of how often each transcript is produced. When all candidates carry the same change, they differ only elsewhere in the gene, and the relative score then reflects those other sites rather than the certainty of the change. The gap between the sample's rank 1 and rank 2, shown below the results table, indicates how clearly the top candidate leads.":
        'Her transkriptin skoru, splice bölgelerinin SpliceAI olasılıklarına göre log-olabilirliğidir: transkriptin kullandığı her aday bölge log(p), kullanmadığı her aday bölge log(1 - p) katkısı yapar. Bu nedenle güçlü bir ekzonun atlanması veya bir intronun tutulması, kullanılmayan güçlü bölgeler oranında skoru düşürür. Göreli skor, gösterilen transkriptlerin skorlarını toplamı 1 olan paylara dönüştürür; bir sıralama ölçütüdür, transkriptlerin hangi sıklıkta üretildiğini göstermez. Tüm adaylar aynı değişikliği taşıdığında yalnızca genin diğer bölgelerinde farklılaşırlar; bu durumda göreli skor, değişikliğin kesinliğinden çok bu diğer bölgeleri yansıtır. Sonuç tablosunun altında gösterilen, örnek dizinin 1. ve 2. sıradaki adayları arasındaki fark, en üstteki adayın ne kadar belirgin biçimde öne çıktığını gösterir.',
    'Evidence vs WT and weak evidence':
        "WT'ye göre kanıt ve zayıf kanıt",
    'For each sample candidate, the largest change in SpliceAI\'s probability (sample minus reference at the same base) at the splice sites where it differs from the WT prediction: sites it uses that the WT prediction does not, and sites it no longer uses. A candidate is marked "weak evidence" when every differing site changed by less than the weak-evidence setting; its difference from the WT prediction then arises from small shifts in probability, often at sites the sequence change barely affects, rather than from a strongly gained or lost site.':
        'Her örnek adayı için, WT tahmininden farklı olduğu splice bölgelerinde (WT tahmininde kullanılmayıp adayda kullanılan ve adayda artık kullanılmayan bölgeler) SpliceAI olasılığındaki en büyük değişimdir (aynı bazda örnek eksi referans). Farklı olduğu tüm bölgelerdeki değişim zayıf kanıt eşiğinin altında kalan adaylar "zayıf kanıt" olarak işaretlenir; bu durumda WT tahmininden farkı, güçlü biçimde kazanılan veya kaybedilen bir bölgeden değil, çoğunlukla dizi değişikliğinden pek etkilenmeyen bölgelerdeki küçük olasılık kaymalarından kaynaklanır.',
    'Gene WT concordance':
        'Gen WT uyumu',
    "How many of this gene's canonical exons the WT prediction reproduces exactly, from a genome-wide run with the default parameters; it is shown on the gene line once the gene is selected. For genes with low concordance, the method's predictions are less reliable.":
        'Varsayılan parametrelerle tüm genom için yapılan bir hesaplamaya göre, WT tahmininin bu genin kanonik ekzonlarından kaçını birebir ürettiğini gösterir; gen seçildiğinde gen satırında görüntülenir. Uyumu düşük genlerde yöntemin tahminleri daha az güvenilirdir.',
    'Local concordance':
        'Yerel uyum',
    'For each change, whether the WT prediction matches the reference in the surrounding exons: the exons on both sides of the affected intron, one further exon on each side, and the introns between them. The exon schematic shades every exon and intron where they differ.':
        'Her değişiklik için, WT tahmininin çevredeki ekzonlarda referansla örtüşüp örtüşmediğini gösterir: etkilenen intronun iki yanındaki ekzonlar, her iki yanda birer ekzon daha ve bunların arasındaki intronlar. Ekzon şemasında, farklılık gösteren tüm ekzon ve intronlar gölgeli olarak gösterilir.',
    'Sequence-specificity control':
        'Dizi özgüllüğü kontrolü',
    'An optional check that runs only when its button is clicked. For each change, the prediction is repeated with random stand-ins in place of the change: an insertion, duplication or inversion is replaced by a random sequence of the same length and dinucleotide composition at the same position, and a deletion by a deletion of the same size at a random position in the same intron. Each change is compared with the sample without it (the other changes kept; the wild type when it is the only change). The tab reports how often the top transcript changes with the stand-ins and their largest delta within 5,000 bases of the change, alongside the values for the real change. If random sequences alter the transcript about as often and as strongly as the real change, the predicted effect is probably related to the position (for example a weak site nearby that any insertion would shift) rather than to the particular sequence. If the real change stands out, its effect is more likely to be specific to its sequence. This is an empirical comparison based on a small number of random draws; it is not a statistical test and does not give a p-value.':
        'Yalnızca ilgili düğmeye tıklandığında çalışan isteğe bağlı bir kontroldür. Her değişiklik için tahmin, değişikliğin yerine rastgele kontrol dizileri konularak yinelenir: insersiyon, duplikasyon veya inversiyon aynı konumda, aynı uzunlukta ve aynı dinükleotit bileşimine sahip rastgele bir diziyle; delesyon ise aynı intronda rastgele bir konumda aynı büyüklükte bir delesyonla değiştirilir. Her değişiklik, örnek dizinin bu değişiklik uygulanmamış haliyle karşılaştırılır (diğer değişiklikler korunur; tek değişiklik söz konusuysa wild type kullanılır). Sekme, kontrol dizileriyle en üstteki transkriptin ne sıklıkta değiştiğini ve değişikliğin 5.000 baz yakınındaki en büyük delta değerlerini, gerçek değişikliğe ait değerlerle birlikte raporlar. Rastgele diziler transkripti gerçek değişiklikle benzer sıklıkta ve güçte etkiliyorsa, öngörülen etki büyük olasılıkla dizinin kendisinden çok konumla ilişkilidir (örneğin herhangi bir insersiyonun kaydırabileceği yakındaki zayıf bir bölge). Gerçek değişiklik belirgin biçimde ayrışıyorsa, etkisinin diziye özgü olma olasılığı daha yüksektir. Bu, az sayıda rastgele denemeye dayanan ampirik bir karşılaştırmadır; istatistiksel bir test değildir ve p değeri üretmez.',
    'Breakpoint sensitivity':
        'Kırılma noktası duyarlılığı',
    'An optional check after Run: each breakpoint of each change is shifted by 10 and 50 bases in both directions, the prediction is repeated, and the result shows whether the top transcript remains the same. A prediction that changes with such small shifts depends on base-exact breakpoints.':
        "Çalıştır'dan sonra kullanılabilen isteğe bağlı bir kontroldür: her değişikliğin kırılma noktaları her iki yöne 10 ve 50 baz kaydırılarak tahmin yinelenir ve en üstteki transkriptin değişip değişmediği raporlanır. Bu kadar küçük kaymalarla değişen bir tahmin, kırılma noktalarının baz düzeyinde doğru olmasına bağlıdır.",
    'RNA check (RT-PCR)':
        'RNA kontrolü (RT-PCR)',
    'For each sample candidate that differs from the WT prediction: the nearest exons on both sides of the change that are unchanged in both the reference transcript and the candidate, where RT-PCR primers could be placed (for example "exon 7 -> exon 9"), and the expected product sizes for the reference and the predicted transcript. Sizes are measured on the cDNAs from the start of the first of these exons to the end of the second; primers placed within the exons shorten both products by the same amount, so the difference between them is preserved. The program does not design primers.':
        'WT tahmininden farklı olan her örnek adayı için: değişikliğin iki yanında, hem referans transkriptte hem de adayda değişmeden korunan ve RT-PCR primerlerinin yerleştirilebileceği en yakın ekzonlar (örneğin "exon 7 -> exon 9") ile referans transkript ve tahmin edilen transkript için beklenen ürün uzunlukları. Uzunluklar cDNA üzerinde, bu ekzonlardan ilkinin başından ikincisinin sonuna kadar ölçülür; primerlerin ekzonların içine yerleştirilmesi her iki ürünü de aynı miktarda kısaltacağından aradaki fark değişmez. Program primer tasarımı yapmaz.',
    'These are predictions for research and hypothesis generation, not measurements. A predicted cDNA may differ from the actual one, and any change must be confirmed with RNA (cDNA) from the sample before a clinical conclusion is drawn.':
        'Bu sonuçlar ölçüm değil, araştırma ve hipotez oluşturma amaçlı tahminlerdir. Tahmin edilen bir cDNA gerçekte oluşandan farklı olabilir; klinik bir sonuca varılmadan önce her değişikliğin örneğe ait RNA (cDNA) ile doğrulanması gerekir.',
    'SpliceAI was trained and evaluated mainly on substitutions and small insertions and deletions. Mid-size insertions and deletions, inversions and duplications fall outside its typical use, so predictions for them are exploratory.':
        'SpliceAI ağırlıklı olarak tek nükleotit değişiklikleri ile küçük insersiyon ve delesyonlar üzerinde eğitilmiş ve değerlendirilmiştir. Orta büyüklükteki insersiyon ve delesyonlar, inversiyonlar ve duplikasyonlar modelin tipik kullanım alanı dışında kaldığından bunlara ilişkin tahminler keşif niteliğindedir.',
    'CNV and SV calls from short-read sequencing often have imprecise breakpoints, whereas SpliceAI requires base-exact breakpoints: a few bases can create or abolish a splice site. Where possible, use breakpoints confirmed by long-read or Sanger sequencing, or run the breakpoint sensitivity check.':
        'Kısa okuma dizilemesinden elde edilen CNV ve SV çağrılarında kırılma noktaları çoğu zaman kesin değildir; SpliceAI ise baz düzeyinde doğru kırılma noktaları gerektirir, çünkü birkaç baz bile bir splice bölgesi oluşturabilir veya ortadan kaldırabilir. Mümkün olduğunda uzun okuma veya Sanger dizilemesiyle doğrulanmış kırılma noktalarının kullanılması ya da kırılma noktası duyarlılığı kontrolünün çalıştırılması önerilir.',
    'SpliceAI can miss some documented deep intronic variants (for example CFTR c.3718-2477C>T, whose 84 bp cryptic exon it does not predict); the absence of a predicted event therefore does not exclude an effect.':
        "SpliceAI, literatürde tanımlanmış bazı derin intronik varyantları öngöremeyebilir (örneğin 84 bp'lik kriptik ekzonunu öngörmediği CFTR c.3718-2477C>T); bu nedenle tahmin edilen bir olayın bulunmaması bir etkinin varlığını dışlamaz.",
    'SpliceAI does not model tissue specificity, the choice of transcription start or polyadenylation site, regulatory factors or isoform usage. Many genes also undergo physiological alternative splicing, which the tab does not label.':
        'SpliceAI doku özgüllüğünü, transkripsiyon başlangıç veya poliadenilasyon bölgesi seçimini, düzenleyici faktörleri ve izoform kullanımını modellemez. Ayrıca birçok gende fizyolojik alternatif splicing bulunur ve sekme bunu ayrıca belirtmez.',
    'The tab shows a notice only when it concerns the current input and affects how the result should be read. A "?" icon lists cautions (for example a sequence that does not cover the whole gene, or a WT prediction that differs from the reference near a change); hover over or click the icon to read them. A red "!" marks a critical problem (for example a sequence that does not align, a reference base that does not match, or a change affecting the ATG or the transcript\'s start or end); only these can prevent Run. In the results, the same icons appear next to the rows concerned. Information that always applies is stated once, here and in the introduction, rather than repeated as a warning each time.':
        'Sekme, bir uyarıyı yalnızca mevcut girdiyle ilgili olduğunda ve sonucun yorumlanmasını etkilediğinde gösterir. "?" simgesi dikkat edilmesi gereken noktaları listeler (örneğin genin tamamını kapsamayan bir dizi veya bir değişikliğin yakınında referanstan farklı olan bir WT tahmini); içeriği görmek için simgenin üzerine gelinebilir veya simgeye tıklanabilir. Kırmızı "!" kritik bir sorunu belirtir (örneğin hizalanamayan bir dizi, referansla eşleşmeyen bir baz veya ATG\'yi ya da transkriptin başlangıç veya bitişini etkileyen bir değişiklik); Çalıştır\'ı yalnızca bu tür sorunlar engelleyebilir. Sonuçlarda aynı simgeler ilgili satırların yanında yer alır. Her durumda geçerli olan bilgiler her seferinde uyarı olarak tekrarlanmaz; bunlar burada ve tanıtımda bir kez belirtilir.',
    "A change in the sample is a stronger hypothesis when four conditions hold together: the WT prediction matches the reference around the change (local concordance); the delta at the changed sites is at least 0.5 (it is not marked as weak evidence); the sample's rank 1 clearly leads rank 2; and the sequence-specificity control shows that random sequences at the same position do not produce the same change.\n\nIt is a weaker hypothesis when the WT prediction is discordant near the change (for example a cryptic exon in that intron present only in the WT prediction), when it is marked as weak evidence, or when random sequences produce the same change.\n\nIn all cases, confirm the change with RNA (cDNA) from the sample; the RNA check column suggests where RT-PCR primers could be placed and which product sizes to expect.":
        "Örnekteki bir değişiklik, aşağıdaki dört koşul birlikte sağlandığında daha güçlü bir hipotez olarak değerlendirilebilir: değişikliğin çevresinde WT tahmini referansla örtüşmektedir (yerel uyum); değişen bölgelerdeki delta en az 0,5'tir (zayıf kanıt olarak işaretlenmemiştir); örnek dizinin 1. sıradaki adayı 2. sıradaki adaydan belirgin biçimde öndedir; dizi özgüllüğü kontrolü, aynı konumdaki rastgele dizilerin aynı değişikliğe yol açmadığını göstermektedir.\n\nWT tahmini değişikliğin yakınında referansla uyumsuz olduğunda (örneğin o intronda yalnızca WT tahmininde bulunan bir kriptik ekzon), aday zayıf kanıt olarak işaretlendiğinde veya rastgele diziler aynı değişikliği oluşturduğunda hipotez daha zayıftır.\n\nHer durumda değişikliğin örneğe ait RNA (cDNA) ile doğrulanması gerekir; RNA kontrolü sütunu, RT-PCR primerlerinin yerleştirilebileceği bölgeleri ve beklenen ürün uzunluklarını önerir.",
    'Ideally, the sequence covers the whole gene. Shorter sequences can also be used: regions not covered are taken from the reference and treated as normal, and every position always receives 5,000 bases of reference context on each side.':
        'Dizinin genin tamamını kapsaması idealdir. Daha kısa diziler de kullanılabilir: kapsanmayan bölgeler referanstan alınır ve normal kabul edilir; her konum için her iki yanda 5.000 bazlık referans bağlamı kullanılır.',
    'default':
        'varsayılan',
    'Always keep in mind:':
        'Her zaman aklınızda tutun:',
    '"Check compatibility" places the sequence on the reference around the selected gene (or within the coordinate window, if one was given) before any prediction is made.\n\nThe summary shows the identity, the part of the reference covered by the sequence, the part of the transcript that is missing (taken from the reference and treated as normal), the orientation in which the sequence was given, and the number of differences. A sequence given in transcript orientation, as exported by Ensembl and NCBI, is the usual case for a minus-strand gene and requires no action. Details contains the technical information: the alignment method, the aligned columns, the positions in the sequence, and where SpliceAI\'s gene interval differs from the canonical transcript. A notice is shown when the sequence clearly fits the other genome build better.\n\nAll differences are listed: substitutions, insertions, deletions, duplications (insertions repeating the adjacent sequence), inversions and complex rearrangements. Insertions and deletions are shifted to their leftmost equivalent position, and each difference is labelled by type. Location gives its position in the canonical transcript, Inserted the inserted bases (shortened in the middle; hover to see all of them, right-click to copy them), and Origin the source of those bases. Only the differences selected under Apply are used to build the sample\'s sequence: by default, those of 10 bp or more and all inversions and duplications, because smaller differences often represent sequencing noise or known polymorphisms. The selection can be changed as needed.\n\n':
        '"Uyumluluğu kontrol et", herhangi bir tahmin yapılmadan önce diziyi seçilen genin çevresinde (koordinat aralığı girilmişse bu aralıkta) referans üzerine yerleştirir.\n\nÖzet; benzerlik oranını, dizinin referansta kapsadığı bölgeyi, transkriptin eksik kalan kısmını (bu kısım referanstan alınır ve normal kabul edilir), dizinin hangi yönde verildiğini ve fark sayısını gösterir. Ensembl ve NCBI\'nin dışa aktardığı biçimde transkript yönünde verilen bir dizi, eksi zincirdeki genler için olağan bir durumdur ve herhangi bir işlem gerektirmez. Ayrıntılar bölümü teknik bilgileri içerir: hizalama yöntemi, hizalanan sütunlar, dizideki konumlar ve SpliceAI gen aralığının kanonik transkriptten farklı olduğu yerler. Dizi diğer genom build\'ine belirgin biçimde daha iyi uyuyorsa bir uyarı gösterilir.\n\nTüm farklar listelenir: tek nükleotit değişiklikleri, insersiyonlar, delesyonlar, duplikasyonlar (komşu diziyi tekrarlayan insersiyonlar), inversiyonlar ve karmaşık yeniden düzenlenmeler. İnsersiyon ve delesyonlar en soldaki eşdeğer konuma kaydırılır ve her fark türüne göre etiketlenir. Konum sütunu farkın kanonik transkriptteki yerini, Eklenen sütunu eklenen bazları (ortası kısaltılarak; tamamını görmek için üzerine gelinebilir, kopyalamak için sağ tıklanabilir), Kaynak sütunu ise bu bazların nereden geldiğini gösterir. Örneğin dizisi yalnızca Uygula sütununda seçili farklarla oluşturulur: varsayılan olarak 10 bp ve üzerindeki farklar ile tüm inversiyon ve duplikasyonlar seçilidir, çünkü daha küçük farklar çoğunlukla dizileme gürültüsü veya bilinen polimorfizmlerdir. Seçim gerektiğinde değiştirilebilir.\n\n',
    'Instead of pasting a sequence, the events themselves can be entered, one per line, and the tab builds the sample\'s sequence from the reference. Supported events: a deletion (chr17:17228330-17228462, optionally followed by the deleted bases, which are then checked), a tandem duplication (a copy immediately after the stretch; write "inverted tandem duplication" or add "inverted" for a reverse-complemented copy), an inversion, or an insertion (a position and the inserted bases, inserted after that position, as in VCF). HGVS g. notation on the chromosome\'s RefSeq accession is also accepted, e.g. NC_000017.11:g.17228330_17228462del, ...dup, ...inv, g.100_101insACGT, g.100_105delinsACG or g.100A>G; the accession must belong to the selected genome build.\n\nCoordinates are 1-based and inclusive, as in VCF files and genome browsers. "Generate sequence" builds the sample\'s sequence for the gene region from the reference and these events and places it in the paste box ("Save as FASTA..." saves it). As the events are known exactly, they go directly into the differences table; the reference is still checked: every event must lie within the gene region, and any reference bases given must match the reference genome.\n\n':
        'Dizi yapıştırmak yerine olaylar her satıra bir tane olacak şekilde girilebilir; sekme örneğin dizisini referanstan oluşturur. Desteklenen olaylar: delesyon (chr17:17228330-17228462; isteğe bağlı olarak ardından silinen bazlar yazılabilir, bu bazlar kontrol edilir), tandem duplikasyon (bölgenin hemen ardından bir kopya; ters-tamamlayıcı kopya için "inverted tandem duplication" yazılabilir veya "inverted" eklenebilir), inversiyon ve insersiyon (bir konum ve eklenen bazlar; VCF\'te olduğu gibi bu konumdan sonra eklenir). Kromozomun RefSeq erişim numarası üzerinden HGVS g. gösterimi de kabul edilir; örneğin NC_000017.11:g.17228330_17228462del, ...dup, ...inv, g.100_101insACGT, g.100_105delinsACG veya g.100A>G. Erişim numarası seçilen genom build\'ine ait olmalıdır.\n\nKoordinatlar, VCF dosyalarında ve genom tarayıcılarında olduğu gibi 1\'den başlar ve uç noktaları kapsar. "Dizi oluştur", gen bölgesine ait örnek dizisini referanstan ve bu olaylardan oluşturarak yapıştırma alanına yerleştirir ("FASTA olarak kaydet..." ile kaydedilebilir). Olaylar kesin olarak bilindiğinden doğrudan farklar tablosuna aktarılır; referans yine de kontrol edilir: her olay gen bölgesi içinde yer almalı ve girilen referans bazları referans genomla eşleşmelidir.\n\n',
    'Optional. A range, such as chr7:117,480,025-117,668,665, restricts the alignment to that part of the reference, which is faster for a short sequence. A single position, such as chr7:117,559,590, makes the check search for the sequence around it: first within twice the sequence length plus 10,000 bases on each side, then across the whole gene if it is not found there. The range found is then entered here.\n\nCoordinates are 1-based and inclusive, as in VCF files and genome browsers. If no gene has been entered, entering coordinates suggests the gene or genes located there.':
        "İsteğe bağlıdır. chr7:117,480,025-117,668,665 gibi bir aralık, hizalamayı referansın bu bölgesiyle sınırlar; bu, kısa diziler için daha hızlıdır. chr7:117,559,590 gibi tek bir konum girildiğinde dizi bu konumun çevresinde aranır: önce her iki yanda dizi uzunluğunun iki katı artı 10.000 baz içinde, burada bulunamazsa genin tamamında. Bulunan aralık daha sonra bu alana yazılır.\n\nKoordinatlar, VCF dosyalarında ve genom tarayıcılarında olduğu gibi 1'den başlar ve uç noktaları kapsar. Gen girilmemişse, koordinat girildiğinde o bölgedeki gen veya genler önerilir.",
    'This tab requires two resources that the VCF tab does not: SpliceAI\'s predictions for the reference genome, stored in a reference score database (about 1.1 GB per genome build), and the reference genome itself. No database is installed on this computer yet; it is provided on the program\'s USB drive and is not distributed through GitHub.\n\n1. The reference genome (hg19.fa, hg38.fa): VCF tab > Settings shows whether it is present; if not, Advanced settings > Download... obtains it from UCSC (about 1 GB per build, once).\n\n2. The database, in one of the following ways:\n• With the USB drive: click "From a folder (USB drive)..." in this tab and select spliceai-refdb\\hg19, then spliceai-refdb\\hg38 on the drive. Every file is checked before use. (Running Setup from the drive does the same.)\n• Without the USB drive: "Compute it on this computer..." computes it with SpliceAI. This takes many hours per genome build and runs in the background; the program remains usable, and the computation resumes after a restart.\n• Without a database: with "Continue without it", each gene is computed when it is analysed (a minute or two, kept for later use).\n\nPurpose of the database: each analysis compares the sample\'s sequence with the reference, and predicting the reference for a whole gene is the slowest step. Once this has been done, each analysis only needs to predict the sample.':
        'Bu sekme, VCF sekmesinde gerekmeyen iki kaynağa ihtiyaç duyar: referans genom için SpliceAI tahminlerinin saklandığı referans skor veritabanı (genom build\'i başına yaklaşık 1,1 GB) ve referans genomun kendisi. Bu bilgisayarda henüz bir veritabanı bulunmamaktadır; veritabanı programın USB belleğiyle birlikte sağlanır, GitHub üzerinden dağıtılmaz.\n\n1. Referans genom (hg19.fa, hg38.fa): VCF sekmesi > Ayarlar bölümünde mevcut olup olmadığı görülebilir; mevcut değilse Gelişmiş ayarlar > İndir... ile UCSC\'den indirilebilir (build başına yaklaşık 1 GB, bir kez).\n\n2. Veritabanı, aşağıdaki yollardan biriyle:\n• USB bellek ile: bu sekmede "Klasörden (USB bellek)..." seçeneğine tıklanarak bellekteki önce spliceai-refdb\\hg19, ardından spliceai-refdb\\hg38 klasörü seçilir. Her dosya kullanılmadan önce kontrol edilir. (Setup\'ın bellekten çalıştırılması da aynı işlemi yapar.)\n• USB bellek olmadan: "Bu bilgisayarda hesapla..." veritabanını SpliceAI ile hesaplar. İşlem genom build\'i başına saatler sürer ve arka planda çalışır; bu sırada program kullanılabilir ve yeniden başlatmanın ardından işlem kaldığı yerden devam eder.\n• Veritabanı olmadan: "Onsuz devam et" seçildiğinde her gen analiz sırasında hesaplanır (bir iki dakika sürer ve sonraki kullanımlar için saklanır).\n\nVeritabanının amacı: her analizde örneğin dizisi referansla karşılaştırılır ve bir genin tamamı için referans tahmini işlemin en uzun süren adımıdır. Bu adım bir kez yapıldığında, sonraki analizlerde yalnızca örnek için tahmin yapılır.',
    'The Sequence-to-cDNA Tab':
        'Sequence to cDNA Sekmesi',
    'Before the First Analysis':
        'İlk Analizden Önce',
    'What It Does':
        'Sekmenin İşlevi',
    'The Parameters':
        'Parametreler',
    'The Results':
        'Sonuçlar',
    'Always Keep in Mind':
        'Göz Önünde Bulundurulması Gerekenler',
    'Notices':
        'Uyarılar',
    'How to Read the Results':
        'Sonuçların Yorumlanması',
    # --- cdna_tab.py
    '{n:,} bp in all: right-click to copy them.':
        'toplam {n:,} bp: kopyalamak için sağ tıklayın.',
    'Right-click to copy.':
        'Kopyalamak için sağ tıklayın.',
    "Don't show again":
        'Bir daha gösterme',
    'Close':
        'Kapat',
    'The introduction to this tab: what it does, every parameter and result, and how to read the result.':
        'Bu sekmenin tanıtımı: ne yaptığı, her parametre ve sonuç, ve sonucun nasıl okunacağı.',
    'Continue without it (genes computed live, a minute or two each)':
        'Onsuz devam et (genler live hesaplanır, her biri bir iki dakika)',
    'Sample sequence':
        'Örnek dizisi',
    'Genome build:':
        "Genom build'i:",
    'choose':
        'seçin',
    'Gene:':
        'Gen:',
    'Gene symbol, e.g. BRCA2':
        'Gen sembolü, ör. BRCA2',
    "SpliceAI's gene symbol or RefSeq's: both are accepted (SpliceAI's annotation uses some older symbols, e.g. FAM132A for C1QTNF12).":
        "SpliceAI'nin ya da RefSeq'in gen sembolü: ikisi de kabul edilir (SpliceAI'nin anotasyonu bazı eski sembolleri kullanır, ör. C1QTNF12 için FAM132A).",
    'Coordinates (optional):':
        'Koordinatlar (isteğe bağlı):',
    'chr7:117,480,025-117,668,665 or chr7:117,559,590':
        'chr7:117,480,025-117,668,665 ya da chr7:117,559,590',
    'A range restricts the alignment to it; a single position is searched around. 1-based, inclusive.':
        "Bir aralık hizalamayı o bölgeyle sınırlar; tek bir konumun çevresinde arama yapılır. 1'den başlar, uçlar dahil.",
    'Coordinates':
        'Koordinatlar',
    'Build this gene now':
        'Bu geni şimdi hesapla',
    'Adds this gene to the reference database first (in the background); otherwise Run computes it live.':
        'Bu geni önce (arka planda) referans veritabanına ekler; aksi halde Çalıştır onu live hesaplar.',
    'Reference score database':
        'Referans skor veritabanı',
    'Paste the sample\'s sequence (raw or FASTA), use "Load FASTA...", or drag a .fa / .fasta / .txt file here.':
        'Örneğin dizisini yapıştırın (düz ya da FASTA), "FASTA yükle..."yi kullanın ya da buraya bir .fa / .fasta / .txt dosyası sürükleyin.',
    'Load FASTA...':
        'FASTA yükle...',
    'Clear':
        'Temizle',
    'Paste sequence':
        'Dizi yapıştır',
    'From coordinates':
        'Koordinatlardan',
    'One event per line. Coordinates are 1-based and inclusive, as in VCF files and genome browsers.':
        "Her satıra bir olay. Koordinatlar VCF dosyalarında ve genom tarayıcılarında olduğu gibi 1'den başlar ve uçları dahildir.",
    'deletion chr17:17228330-17228462\ntandem duplication chr17:17228330-17228462 (add "inverted" for an inverted one)\ninversion chr17:17228330-17228462\ninsertion chr17:17228330 ACGTACGT (inserted after that position)\nNC_000017.11:g.17228330_17228462del':
        'deletion chr17:17228330-17228462\ntandem duplication chr17:17228330-17228462 (inverte için "inverted" ekleyin)\ninversion chr17:17228330-17228462\ninsertion chr17:17228330 ACGTACGT (o konumdan sonra eklenir)\nNC_000017.11:g.17228330_17228462del',
    'Generate sequence':
        'Dizi oluştur',
    "The sample sequence of the gene's region: the reference with these events. It goes into the paste box, and the events into the differences table.":
        'Genin bölgesinin örnek dizisi: bu olayların uygulandığı referans. Yapıştırma kutusuna, olaylar da farklar tablosuna gider.',
    'Save as FASTA...':
        'FASTA olarak kaydet...',
    'Compatibility':
        'Uyumluluk',
    'Check compatibility':
        'Uyumluluğu kontrol et',
    'Cautions about this input':
        'Bu girdiyle ilgili dikkat notları',
    'Critical problems':
        'Kritik sorunlar',
    'Details':
        'Ayrıntılar',
    "Run anyway (the change touches the ATG or the transcript's start or end; the predicted transcripts may lose that anchor)":
        "Yine de çalıştır (değişiklik ATG'ye ya da transkriptin başına veya sonuna dokunuyor; tahmin edilen transkriptler bu çapayı kaybedebilir)",
    'Parameters':
        'Parametreler',
    'Site threshold:':
        'Splice bölgesi eşiği:',
    'Minimum internal exon:':
        'En kısa iç ekzon:',
    ' nt':
        ' nt',
    'Minimum intron:':
        'En kısa intron:',
    'Apply by default from:':
        'Varsayılan uygulama sınırı:',
    ' bp':
        ' bp',
    'Differences this long or longer (and every inversion and duplication) are ticked under Apply when the check finds them.':
        'Bu uzunlukta ya da daha uzun farklar (ve her inversiyon ve duplikasyon) kontrol onları bulduğunda Uygula altında işaretlenir.',
    'Advanced':
        'Gelişmiş',
    'Transcripts per sequence:':
        'Dizi başına transkript:',
    'Weak evidence below delta:':
        'Zayıf kanıt delta sınırı:',
    'A sample candidate whose differing splice sites all changed by less than this is marked "weak evidence" (see Reading the results). It doesn\'t change the ranking.':
        'Farklı splice bölgelerinin hepsi bundan daha az değişen bir örnek adayı "zayıf kanıt" olarak işaretlenir (bkz. Sonuçları okuma). Sıralamayı değiştirmez.',
    'Random stand-ins:':
        'Rastgele kontrol dizileri:',
    'How many random stand-ins per change the sequence-specificity control uses.':
        'Dizi özgüllüğü kontrolünün her değişiklik için kullandığı rastgele kontrol dizisi sayısı.',
    'Run':
        'Çalıştır',
    'Cancel':
        'İptal',
    'Results':
        'Sonuçlar',
    'WT alternatives (2 to 5)':
        'WT alternatifleri (2-5)',
    "Also shows the wild type's alternative transcripts, ranks 2 to 5.":
        "Wild type'ın 2-5. sıradaki alternatif transkriptlerini de gösterir.",
    'Export CSV...':
        'CSV dışa aktar...',
    'The rows shown, as a CSV or tab-separated file.':
        'Gösterilen satırlar, CSV ya da sekmeyle ayrılmış dosya olarak.',
    'Download results...':
        'Sonuçları indir...',
    'Everything, in one file this tab (and the viewer) can open again without SpliceAI, the FASTA or the database.':
        'Her şey, bu sekmenin (ve görüntüleyicinin) SpliceAI, FASTA ya da veritabanı olmadan yeniden açabileceği tek bir dosyada.',
    'Open calculated results...':
        'Hesaplanmış sonuçları aç...',
    'Opens a results file saved with "Download results" -- or drag it onto this tab.':
        '"Sonuçları indir" ile kaydedilmiş bir sonuç dosyasını açar -- ya da dosyayı bu sekmeye sürükleyin.',
    'Open viewer':
        'Görüntüleyiciyi aç',
    'The base-level view, the exon schematic and the sequence comparison -- or click a row.':
        'Baz düzeyinde görünüm, ekzon şeması ve dizi karşılaştırması -- ya da bir satıra tıklayın.',
    'Export bedGraph...':
        'bedGraph dışa aktar...',
    'Reference, sample and delta tracks for IGV or the UCSC browser.':
        'IGV ya da UCSC tarayıcısı için referans, örnek ve delta izleri.',
    'Optional: repeats the prediction with random stand-ins for each change (see Reading the results). Runs only when clicked.':
        'İsteğe bağlı: tahmini her değişiklik için rastgele kontrol dizileriyle tekrarlar (bkz. Sonuçları okuma). Yalnızca tıklanınca çalışır.',
    'Optional: moves each breakpoint by 10 and 50 bases either way and says whether the top transcript changes. Runs only when clicked.':
        'İsteğe bağlı: her kırılma noktasını iki yöne 10 ve 50 baz kaydırır ve en üstteki transkriptin değişip değişmediğini söyler. Yalnızca tıklanınca çalışır.',
    "Why the VCF tab's scores differ":
        'VCF sekmesinin skorları neden farklı',
    'Sequence to cDNA':
        'Sequence to cDNA',
    'Choose the genome build and the gene first.':
        "Önce genom build'ini ve geni seçin.",
    'Building the sample sequence...':
        'Örnek dizisi oluşturuluyor...',
    'Generated from {n} event(s) given as coordinates ({bases:,} bases)':
        'Koordinat olarak verilen {n} olaydan oluşturuldu ({bases:,} baz)',
    '{n} event(s) applied to {region}; the sequence is in the paste box.':
        '{n} olay {region} bölgesine uygulandı; dizi yapıştırma kutusunda.',
    'Save as FASTA':
        'FASTA olarak kaydet',
    'FASTA (*.fa *.fasta)':
        'FASTA (*.fa *.fasta)',
    'Looking for genes there...':
        'Oradaki genler aranıyor...',
    "No gene of SpliceAI's {build} gene table at {where}.":
        "SpliceAI'nin {build} gen tablosunda {where} konumunda gen yok.",
    'Genes at {where}: {names}; {gene} is chosen (change it in the list). ':
        '{where} konumundaki genler: {names}; {gene} seçildi (listeden değiştirebilirsiniz). ',
    'Gene at {where}: {names}; {gene} is chosen. ':
        '{where} konumundaki gen: {names}; {gene} seçildi. ',
    'Looking the gene up...':
        'Gen aranıyor...',
    'no NM_':
        'NM_ yok',
    'No gene "{name}" in SpliceAI\'s {build} gene table.':
        'SpliceAI\'nin {build} gen tablosunda "{name}" adında gen yok.',
    '+ strand':
        '+ zincir',
    '− strand':
        '− zincir',
    'no NM_ transcript':
        'NM_ transkripti yok',
    'The locus is shared with {genes}.':
        'Bu lokus {genes} ile ortak.',
    ' and {n} more':
        ' ve {n} tane daha',
    'WT concordance for this gene: {same}/{exons} exons':
        'Bu genin WT uyumu: {same}/{exons} ekzon',
    ', the whole transcript exactly':
        ', transkriptin tamamı birebir',
    '; predicted differently: exons {list}':
        '; farklı tahmin edilen ekzonlar: {list}',
    '; predicted differently: exon {list}':
        '; farklı tahmin edilen ekzon: {list}',
    ' (default parameters)':
        ' (varsayılan parametreler)',
    'Reference values: from the database.':
        'Referans değerleri: veritabanından.',
    'Reference values: computed live (not in the database), about {time} for its {bases:,} bases, once; then cached.':
        'Referans değerleri: live hesaplanır (veritabanında yok), {bases:,} bazı için yaklaşık {time}, bir kez; sonra önbellekte tutulur.',
    'Building it now (Reference database at the top).':
        'Şimdi hesaplanıyor (üstteki Reference DB).',
    "Load the sample's sequence":
        'Örneğin dizisini yükle',
    'Sequence files (*.fa *.fasta *.fna *.txt);;All files (*)':
        'Dizi dosyaları (*.fa *.fasta *.fna *.txt);;Tüm dosyalar (*)',
    'Load FASTA':
        'FASTA yükle',
    "{path} can't be read: {error}":
        '{path} okunamıyor: {error}',
    "Can't use this file":
        'Bu dosya kullanılamıyor',
    "{name} can't be loaded here.\n\nDrag in a sequence (.fa, .fasta, .txt) or a results file (*{ext}).":
        '{name} buraya yüklenemez.\n\nBir dizi (.fa, .fasta, .txt) ya da bir sonuç dosyası (*{ext}) sürükleyin.',
    'Reading and aligning the sequence...':
        'Dizi okunuyor ve hizalanıyor...',
    "✗ The sequence can't be used.":
        '✗ Dizi kullanılamıyor.',
    '✓ Checked.':
        '✓ Kontrol edildi.',
    'Copy inserted sequence ({n:,} bp)':
        'Eklenen diziyi kopyala ({n:,} bp)',
    'Estimated run time: about {time} ({n} stretch(es) near the edits predicted again).':
        'Tahmini süre: yaklaşık {time} (değişikliklerin yakınındaki {n} bölge yeniden tahmin edilir).',
    'Estimated run time: about {time} ({n} stretch(es) near the edits predicted again; the reference computed live first).':
        'Tahmini süre: yaklaşık {time} (değişikliklerin yakınındaki {n} bölge yeniden tahmin edilir; önce referans live hesaplanır).',
    'Run is blocked: {n} critical problem(s) -- the red ! on the right says which.':
        'Çalıştır engellendi: {n} kritik sorun -- sağdaki kırmızı ! hangileri olduğunu söyler.',
    "✗ A change touches the ATG or the transcript's start or end (the red !): tick Run anyway to go on.":
        "✗ Bir değişiklik ATG'ye ya da transkriptin başına veya sonuna dokunuyor (kırmızı !): devam etmek için Yine de çalıştır'ı işaretleyin.",
    'Checked.':
        'Kontrol edildi.',
    '{n} caution(s): the ? on the right.':
        '{n} dikkat notu: sağdaki ?.',
    'Starting':
        'Başlıyor',
    '{step} done in {time}.':
        '{step}: {time} içinde bitti.',
    'Sequence-specificity control (an empirical comparison, not a p-value): ':
        'Dizi özgüllüğü kontrolü (deneysel bir karşılaştırma, p değeri değil): ',
    'Breakpoint sensitivity (each breakpoint moved by 10 and 50 bases either way): ':
        'Kırılma noktası duyarlılığı (her kırılma noktası iki yöne 10 ve 50 baz kaydırıldı): ',
    'Done in {time}.':
        '{time} içinde bitti.',
    "The reference values were computed live (this gene isn't in the reference score database yet); the results are the same, only slower.":
        'Referans değerleri live hesaplandı (bu gen henüz referans skor veritabanında bulunmuyor); sonuçlar aynıdır, yalnızca işlem daha uzun sürer.',
    'Failed.':
        'Başarısız.',
    'Run failed':
        'Çalıştırma başarısız',
    'Cancelled.':
        'İptal edildi.',
    'The {n} variant(s) still queued are not run.':
        'Kuyrukta kalan {n} varyant çalıştırılmayacak.',
    '{n} variant(s), results in {folder}:':
        '{n} varyant, sonuçlar {folder} klasöründe:',
    'From the VCF tab: {variant} in {gene}':
        'VCF sekmesinden: {gene} içinde {variant}',
    '({i} of {n})':
        '({i}/{n})',
    'Analyse transcript effect':
        'Transkript etkisini analiz et',
    "{variant} isn't run: {why}\n\nThe check's notices are next to Check compatibility; once they allow it, Run starts it.":
        "{variant} çalıştırılmadı: {why}\n\nKontrolün uyarıları Uyumluluğu kontrol et'in yanında; izin verdiklerinde Çalıştır onu başlatır.",
    'The variant':
        'Varyant',
    'Applied together with {n} other variant(s) of the VCF on one haplotype (their phase is unknown): ':
        "VCF'teki {n} başka varyantla tek bir haplotip üzerinde birlikte uygulandı (fazları bilinmiyor): ",
    'within the rounding tolerance (±0.005)':
        'yuvarlama toleransı (±0,005) içinde',
    'differs beyond the rounding tolerance (the ? says why): ':
        'yuvarlama toleransının dışında farklı (? nedenini söyler): ',
    '<b>From the VCF tab</b>: {variant} in {gene}; genotype {genotype}':
        '<b>VCF sekmesinden</b>: {gene} içinde {variant}; genotip {genotype}',
    'Cancelling...':
        'İptal ediliyor...',
    ', about {time} left':
        ', yaklaşık {time} kaldı',
    'Elapsed {time}':
        'Geçen süre {time}',
    'partly built':
        'kısmen hesaplandı',
    'being built':
        'hesaplanıyor',
    'missing':
        'yok',
    'not usable':
        'kullanılamaz',
    'checking':
        'kontrol ediliyor',
    'Database {build}: {state}':
        'Veritabanı {build}: {state}',
    '(version {version})':
        '(sürüm {version})',
    'cache {mb:,.1f} MB':
        'önbellek {mb:,.1f} MB',
    '<b>Applied events</b> (click one to see it in the viewer):':
        '<b>Uygulanan olaylar</b> (görüntüleyicide görmek için birine tıklayın):',
    'Export CSV':
        'CSV dışa aktar',
    'CSV (*.csv);;Tab-separated (*.tsv)':
        'CSV (*.csv);;Sekmeyle ayrılmış (*.tsv)',
    'Download results':
        'Sonuçları indir',
    'Sequence-to-cDNA results (*{ext})':
        'Sequence to cDNA sonuçları (*{ext})',
    'Open calculated results':
        'Hesaplanmış sonuçları aç',
    'Opened {name}: {gene} on {build}, made {created}.':
        '{name} açıldı: {gene}, {build}, oluşturulma {created}.',
    'Folder for the bedGraph tracks':
        'bedGraph izleri için klasör',
    # --- cdna_worker.py
    "SpliceAI is not installed, so the sample's sequence cannot be predicted. To install it, click the SpliceAI button at the top of the window.":
        'SpliceAI kurulu olmadığından örneğin dizisi için tahmin yapılamıyor. Kurulum için pencerenin üstündeki SpliceAI düğmesine tıklayın.',
    # --- coord_events.py
    'Line {n} ({line}): not an event like "deletion chr17:17228330-17228462", "insertion chr17:17228330 ACGT" or NC_000017.11:g.17228330_17228462del.':
        'Satır {n} ({line}): "deletion chr17:17228330-17228462", "insertion chr17:17228330 ACGT" ya da NC_000017.11:g.17228330_17228462del gibi bir olay değil.',
    'Line {n}: an insertion is a position and the inserted bases, e.g. "insertion chr17:17228330 ACGT" (inserted after that position).':
        'Satır {n}: insersiyon bir konum ve eklenen bazlardır, ör. "insertion chr17:17228330 ACGT" (o konumdan sonra eklenir).',
    'Line {n}: only a deletion takes bases after its coordinates (the deleted ones).':
        'Satır {n}: koordinatlarından sonra baz alan tek olay delesyondur (silinen bazlar).',
    'Line {n}: the start has to be at least 1 and not after the end.':
        'Satır {n}: başlangıç en az 1 olmalı ve bitişten sonra olmamalı.',
    "Line {n}: {acc} isn't a GRCh37 or GRCh38 chromosome accession.":
        'Satır {n}: {acc} bir GRCh37 ya da GRCh38 kromozom erişim numarası değil.',
    'Line {n}: {acc} is a {assembly} chromosome, but {build} is chosen.':
        'Satır {n}: {acc} bir {assembly} kromozomu, ama {build} seçili.',
    'Line {n}: a substitution is at one position.':
        'Satır {n}: tek baz değişikliği tek bir konumdadır.',
    'Line {n}: an HGVS insertion names the two bases it lies between (g.100_101insACGT) and the inserted bases.':
        'Satır {n}: bir HGVS insersiyonu arasında yer aldığı iki bazı (g.100_101insACGT) ve eklenen bazları belirtir.',
    'Line {n}: delins needs the inserted bases.':
        'Satır {n}: delins için eklenen bazlar gerekir.',
    'Line {n}: {bases} bases given for {positions:,} positions.':
        'Satır {n}: {positions:,} konum için {bases} baz verildi.',
    'Line {n}: chr{chrom}, but the gene is on chr{gene_chrom}.':
        'Satır {n}: chr{chrom}, ama gen chr{gene_chrom} üzerinde.',
    "Line {n}: {end:,} is past the chromosome's end ({length:,}).":
        'Satır {n}: {end:,}, kromozomun sonunu ({length:,}) geçiyor.',
    'Line {n}: {start:,}..{end:,} is outside the region analysed for this gene ({r0:,}..{r1:,}).':
        'Satır {n}: {start:,}..{end:,}, bu gen için analiz edilen bölgenin ({r0:,}..{r1:,}) dışında.',
    'Line {n}: the reference has N (unknown bases) in {start:,}..{end:,}.':
        'Satır {n}: referansta {start:,}..{end:,} arasında N (bilinmeyen baz) var.',
    'Line {n}: the reference at {where} is {ref}, not {given}.':
        'Satır {n}: {where} konumunda referans {ref}, {given} değil.',
    'The {kind_a} at {a:,} and the {kind_b} at {b:,} overlap.':
        '{a:,} konumundaki {kind_a} ile {b:,} konumundaki {kind_b} örtüşüyor.',
    # --- download_dialog.py
    'Downloading':
        'İndiriliyor',
    'Decompressing':
        'Açılıyor',
    'Building index':
        'İndeks oluşturuluyor',
    'Download reference FASTA':
        'Referans FASTA indir',
    'Build:':
        'Build:',
    'Save to folder:':
        'Kaydedilecek klasör:',
    'Browse...':
        'Gözat...',
    "Downloads about 1 GB (compressed) from UCSC's goldenPath, decompresses it to about 3 GB and builds the .fai index. Depending on the connection, this may take some time; the progress bar below shows the progress.":
        "UCSC goldenPath'ten yaklaşık 1 GB'lık sıkıştırılmış dosyayı indirir, yaklaşık 3 GB olarak açar ve .fai indeksini oluşturur. Bağlantı hızına bağlı olarak işlem zaman alabilir; ilerleme aşağıdaki çubukta gösterilir.",
    'Download':
        'İndir',
    'Select destination folder':
        'Hedef klasörü seçin',
    'No folder selected':
        'Klasör seçilmedi',
    'Choose a destination folder first.':
        'Önce bir hedef klasör seçin.',
    'Starting...':
        'Başlıyor...',
    'Done: {path}':
        'Bitti: {path}',
    'Download failed':
        'İndirme başarısız',
    # --- exon_view.py
    'WT prediction differs from the reference in this region; interpret the sample prediction here with caution.':
        'Bu bölgede WT tahmini referanstan farklı; buradaki örnek tahminini dikkatle yorumlayın.',
    "Transcript orientation (5'→3'), {strand}: exon 1 on the left":
        "Transkript yönü (5'→3'), {strand}: ekzon 1 solda",
    "Transcript orientation (5'→3'), {strand}":
        "Transkript yönü (5'→3'), {strand}",
    'Genomic orientation, {strand}':
        'Genomik yön, {strand}',
    'Genomic orientation, {strand}: exon 1 on the right':
        'Genomik yön, {strand}: ekzon 1 sağda',
    "SpliceAI (WT) at the reference's splice sites: {sites} (site threshold {threshold:g})":
        'Referansın splice bölgelerinde SpliceAI (WT): {sites} (splice bölgesi eşiği {threshold:g})',
    'WT prediction differs from the reference here: {kind} {k} ({what})':
        'Burada WT tahmini referanstan farklı: {kind} {k} ({what})',
    'no change vs the reference':
        'referansa göre değişiklik yok',
    'no change vs WT predicted 1':
        "WT predicted 1'e göre değişiklik yok",
    'exon {n} (dup): a copy of reference exon {n}':
        'ekzon {n} (dup): referans ekzon {n} kopyası',
    'an exon in inserted sequence':
        'eklenen dizide bir ekzon',
    '{row}: {what}, {n:,} nt, in the inserted bases':
        '{row}: {what}, {n:,} nt, eklenen bazlarda',
    'reference exon {n}':
        'referans ekzon {n}',
    'no reference exon (novel)':
        'referans ekzon yok (yeni)',
    '{row}, exon {n} ({same})':
        '{row}, ekzon {n} ({same})',
    '{row}: exon {n} of {total}':
        '{row}: ekzon {n}/{total}',
    'reference {where} ({n:,} nt)':
        'referans {where} ({n:,} nt)',
    'matches reference exon {n}':
        'referans ekzon {n} ile eşleşiyor',
    'haplotype {start:,}-{end:,}':
        'haplotip {start:,}-{end:,}',
    'vs the reference ({accession}): {status}':
        'referansa göre ({accession}): {status}',
    'vs WT predicted 1: {status}':
        "WT predicted 1'e göre: {status}",
    'change vs the reference: {change}':
        'referansa göre değişiklik: {change}',
    # --- first_launch_dialog.py
    'Welcome to SpliceAI Variant Scoring.':
        "SpliceAI Variant Scoring'e hoş geldiniz.",
    "This program was developed so that users who wish to work with SpliceAI can examine their variants without Linux, WSL or the command line, and so that SpliceAI's raw output can be brought into a form suitable for analysis.":
        'Bu program, SpliceAI kullanmak isteyen kullanıcının varyantları Linux, WSL ya da komut satırı gerektirmeden inceleyebilmesi ve SpliceAI ham verilerinin analize uygun hale getirilebilmesi adına yazıldı.',
    'About SpliceAI':
        'SpliceAI Hakkında',
    'SpliceAI is a deep-learning model developed by Illumina that predicts whether a genetic variant alters normal RNA splicing. This program provides an interface to the model; it does not modify or replace it.':
        'SpliceAI, Illumina tarafından geliştirilmiş ve bir genetik varyantın normal RNA splicing sürecini etkileyip etkilemediğini öngören bir derin öğrenme modelidir. Bu program, söz konusu model için bir arayüz sunar; modelde herhangi bir değişiklik yapmaz ve modelin yerini almaz.',
    'About SnpEff':
        'SnpEff Hakkında',
    'About Sequence to cDNA':
        'Sequence to cDNA Hakkında',
    'The Sequence to cDNA tab takes a sample sequence of one gene, or variants entered as coordinates, and uses SpliceAI to predict how it will be spliced. This makes it possible to compare the reference exons, the normal splicing predicted for the reference sequence, and the predicted effect of the variant.':
        "Sequence to cDNA sekmesi, bir gene ait örnek diziyi ya da koordinatla girilen varyantları alır ve bu dizinin nasıl splice edileceğini SpliceAI ile öngörür. Böylece referans ekzonlar, referans dizi için öngörülen normal splicing ve varyantın öngörülen etkisi karşılaştırılabilir.",
    'Terms of Use':
        'Kullanım Koşulları',
    'This software is provided for academic and research use only, without warranty of any kind, express or implied. It is not a diagnostic tool, and its output must not be used for clinical decision-making without independent verification by qualified personnel. SpliceAI itself is licensed separately by Illumina, Inc.; use of this program does not grant any rights under that license, and the user is solely responsible for compliance with it. The developer assumes no liability for any use of this software or reliance on its output. By using this program, the user accepts all associated risk.':
        'Bu yazılım, yalnızca akademik ve araştırma amaçlı kullanım için ve açık ya da örtük herhangi bir garanti verilmeksizin sunulmaktadır. Yazılım bir tanı aracı değildir; çıktıları, yetkin kişilerce bağımsız olarak doğrulanmadan klinik karar süreçlerinde kullanılmamalıdır. SpliceAI, Illumina, Inc. tarafından ayrıca lisanslanmaktadır; bu programın kullanılması söz konusu lisans kapsamında herhangi bir hak doğurmaz ve lisans koşullarına uyulması kullanıcının sorumluluğundadır. Geliştirici, yazılımın kullanımından veya çıktılarına dayanılarak yapılan işlemlerden doğabilecek sonuçlardan sorumlu tutulamaz. Programı kullanan kişi, kullanımdan doğabilecek tüm riskleri kabul etmiş sayılır.',
    'Bug reports and suggestions':
        'Hata bildirimleri ve öneriler',
    'Welcome to SpliceAI Variant Scoring':
        "SpliceAI Variant Scoring'e hoş geldiniz",
    'Opens your e-mail program':
        'E-posta programınızı açar',
    'Read the full license...':
        'Lisansın tamamını oku...',
    'Decline':
        'Reddet',
    'I have read and accept the terms':
        'Koşulları okudum ve kabul ediyorum',
    # --- gene_import_dialog.py
    'Import gene list':
        'Gen listesi içe aktar',
    'Saved lists:':
        'Kayıtlı listeler:',
    'Load':
        'Yükle',
    'Save...':
        'Kaydet...',
    'Delete':
        'Sil',
    'Paste gene symbols below, one per line or separated by commas or spaces (for example, copied from an HPO phenotype gene list):':
        'Gen sembollerini aşağıya yapıştırın; her satıra bir sembol yazılabilir veya semboller virgül ya da boşlukla ayrılabilir (örneğin bir HPO fenotip gen listesinden kopyalanmış olarak):',
    'BRCA1, BRCA2\nANKRD26\nTP53 ...':
        'BRCA1, BRCA2\nANKRD26\nTP53 ...',
    'Check':
        'Kontrol et',
    'Not recognized (check for typos or non-gene text):':
        'Tanınmayanlar (yazım hatası ya da gen olmayan metin olabilir):',
    'Recognized:':
        'Tanınanlar:',
    'Apply Filter':
        'Filtreyi uygula',
    'No gene symbols detected in the pasted text.':
        'Yapıştırılan metinde gen sembolü bulunamadı.',
    'SpliceAI not installed':
        'SpliceAI kurulu değil',
    "Checking gene symbols requires SpliceAI to be installed (it provides the reference gene list used for the check).\n\nTo install it, click the SpliceAI status button in the main window's top bar.":
        "Gen sembollerinin kontrol edilebilmesi için SpliceAI'nin kurulu olması gerekir (kontrolde kullanılan referans gen listesi SpliceAI ile birlikte gelir).\n\nKurulum için ana pencerenin üst çubuğundaki SpliceAI durum düğmesine tıklayın.",
    'Check failed':
        'Kontrol başarısız',
    '{n} entries detected -- {valid} recognized, {invalid} not recognized.':
        '{n} giriş bulundu -- {valid} tanındı, {invalid} tanınmadı.',
    '(none)':
        '(yok)',
    'Save gene list':
        'Gen listesini kaydet',
    'Name for this list:':
        'Bu listenin adı:',
    'Overwrite saved list?':
        'Kayıtlı listenin üzerine yazılsın mı?',
    'A saved list named "{name}" already exists. Overwrite it?':
        '"{name}" adında kayıtlı bir liste zaten var. Üzerine yazılsın mı?',
    'Saved':
        'Kaydedildi',
    'Saved "{name}" ({n} gene(s)).':
        '"{name}" kaydedildi ({n} gen).',
    'Delete saved list?':
        'Kayıtlı liste silinsin mi?',
    'Delete saved list "{name}"?':
        '"{name}" kayıtlı listesi silinsin mi?',
    # --- license_dialog.py
    'Third-Party Licenses -- SpliceAI Variant Scoring':
        'Üçüncü taraf lisansları -- SpliceAI Variant Scoring',
    '<b>Third-Party Licenses</b><br>The license terms of SpliceAI, SnpEff, NCBI MANE and the UCSC reference genome.':
        '<b>Üçüncü Taraf Lisansları</b><br>SpliceAI, SnpEff, NCBI MANE ve UCSC referans genomuna ait lisans koşulları.',
    # --- main_window.py
    '⏸ Paused -- ':
        '⏸ Duraklatıldı -- ',
    "raw: the model's unfiltered delta scores.\n\nmasked: predictions that do not match an annotated splice boundary in the expected direction are set to zero, to suppress biologically implausible predictions. For the same variant, masked scores are never higher than raw scores.":
        'raw: modelin filtrelenmemiş delta skorlarıdır.\n\nmasked: beklenen yönde anotasyonlu bir splice sınırıyla örtüşmeyen tahminler, biyolojik olarak olası olmayan sonuçları bastırmak amacıyla sıfırlanır. Aynı varyant için masked skorlar hiçbir zaman raw skorlardan yüksek değildir.',
    "Column abbreviations:\n  • DS = Delta Score: the predicted probability (0-1) that the variant creates or disrupts a splice site of the given type\n  • DP = Delta Position: the distance in base pairs between the variant and the predicted splicing change (negative or positive, indicating the direction)\n  • AG = Acceptor Gain: a new splice acceptor site is predicted\n  • AL = Acceptor Loss: the loss of an existing splice acceptor site is predicted\n  • DG = Donor Gain: a new splice donor site is predicted\n  • DL = Donor Loss: the loss of an existing splice donor site is predicted\n\nTogether these form four score/position column pairs: DS_AG / DP_AG, DS_AL / DP_AL, DS_DG / DP_DG and DS_DL / DP_DL.\n\nmax_score is the largest of the four scores and summarises the overall likelihood that the variant affects splicing.\n\nCommonly used interpretation thresholds (conventions, not diagnostic criteria): the original SpliceAI publication (Jaganathan et al., Cell 2019) and many laboratories use about 0.2 as a sensitive cut-off, about 0.5 as a standard threshold and about 0.8 as a high-confidence threshold.\n\nA high score indicates a predicted splicing effect; it does not by itself establish clinical pathogenicity. Variant interpretation also requires the ACMG criteria, clinical correlation and other supporting evidence.\n\nAllele Fraction is read from the VCF's AF field (the per-sample FORMAT value is preferred to the INFO value). If there is no AF, it is calculated from AD as alt / (ref + alt); if neither is present, the column is left blank. It is not needed for scoring.":
        "Sütun kısaltmaları:\n  • DS = Delta Score: varyantın ilgili türde bir splice bölgesi oluşturma veya bozma olasılığı (0-1)\n  • DP = Delta Position: öngörülen splicing değişikliğinin varyanta baz çifti cinsinden uzaklığı (negatif veya pozitif değer yönü gösterir)\n  • AG = Acceptor Gain: yeni bir splice akseptör bölgesi oluşması öngörülür\n  • AL = Acceptor Loss: mevcut bir splice akseptör bölgesinin kaybı öngörülür\n  • DG = Donor Gain: yeni bir splice donör bölgesi oluşması öngörülür\n  • DL = Donor Loss: mevcut bir splice donör bölgesinin kaybı öngörülür\n\nBunlar birlikte dört skor/konum sütun çiftini oluşturur: DS_AG / DP_AG, DS_AL / DP_AL, DS_DG / DP_DG ve DS_DL / DP_DL.\n\nmax_score, dört skorun en büyüğüdür ve varyantın splicing'i etkileme olasılığını tek bir değerle özetler.\n\nYaygın olarak kullanılan yorum eşikleri (tanı ölçütü değil, kabul görmüş uygulamalardır): SpliceAI'nin özgün yayını (Jaganathan ve ark., Cell 2019) ve birçok laboratuvar, duyarlı bir eşik olarak yaklaşık 0,2'yi, standart eşik olarak yaklaşık 0,5'i ve yüksek güvenilirlik eşiği olarak yaklaşık 0,8'i kullanmaktadır.\n\nYüksek bir skor öngörülen bir splicing etkisine işaret eder; tek başına klinik patojenite anlamı taşımaz. Varyant yorumu için ayrıca ACMG ölçütleri, klinik bulgularla değerlendirme ve diğer destekleyici kanıtlar gereklidir.\n\nAllele Fraction, VCF dosyasının AF alanından okunur (örneğe özgü FORMAT değeri, INFO değerine tercih edilir). AF bulunmuyorsa AD alanından alt / (ref + alt) olarak hesaplanır; ikisi de yoksa sütun boş bırakılır. Skorlama için gerekli değildir.",
    'The SpliceAI authors (Jaganathan et al., Cell 2019) have precomputed scores for every possible single-nucleotide variant and for common small indels (1 bp insertions, 1-4 bp deletions). These files are optional and can be downloaded from Illumina BaseSpace (basespace.illumina.com/s/otSPW8hnhaZR) under Illumina\'s terms. Scores are available for hg19 and hg38, in raw and masked mode; for one build, the masked SNV file is about 29 GB and the masked indel file about 69 GB.\n\nWhen the program can read these files, variants already in them are not scored again: their scores are taken from the files, which shortens the run. Without them ("I don\'t have precomputed data") every variant is scored live; this takes longer, but the results are the same.\n\nFile names: spliceai_scores.{raw,masked}.{snv,indel}.hg{19,38}.vcf.gz, each with its .tbi index. Keep the files compressed as downloaded (bgzip, with the .tbi index); do not decompress them. Then select "Use precomputed data" and choose the folder that contains them.':
        'SpliceAI\'ın geliştiricileri (Jaganathan ve ark., Cell 2019), olası tüm tek nükleotit varyantları ve yaygın küçük indel\'ler (1 bp insersiyon, 1-4 bp delesyon) için skorları önceden hesaplamıştır. Bu dosyalar isteğe bağlıdır; Illumina\'nın koşullarına uyularak Illumina BaseSpace\'ten (basespace.illumina.com/s/otSPW8hnhaZR) indirilebilir. hg19 ve hg38 için hem raw hem masked skorlar bulunur; bir build için masked SNV dosyası yaklaşık 29 GB, masked indel dosyası yaklaşık 69 GB\'tır.\n\nProgram bu dosyalara erişebildiğinde, dosyalarda bulunan varyantlar yeniden skorlanmaz; skorları dosyadan alınır ve işlem süresi kısalır. Dosyalar olmadan ("Precomputed verim yok") tüm varyantlar live skorlanır; süre uzar, ancak sonuçlar aynıdır.\n\nDosya adları: spliceai_scores.{raw,masked}.{snv,indel}.hg{19,38}.vcf.gz ve her birinin .tbi indeksi. Dosyalar indirildikleri sıkıştırılmış biçimde (bgzip ve .tbi indeksiyle) tutulmalı, açılmamalıdır. Ardından "Precomputed veri kullan" seçilir ve dosyaların bulunduğu klasör gösterilir.',
    "Scores only the variants within the selected genes, so a whole genome or exome takes minutes instead of hours. A variant is considered to lie within a gene when it falls in SpliceAI's interval for that gene (the same rule that assigns a score to a gene), and only rows for those genes are kept.\n\nRun is enabled once one of the two is chosen: All genes scores every variant; Choose gene list... accepts a pasted list (for example an HPO phenotype list or a gene panel) and can save it under a name for later use. The last list is offered again.\n\nThis differs from the gene-set filter above the results table, which only hides rows of a run in which all variants were scored.":
        "Yalnızca seçilen genlerin içindeki varyantlar skorlanır; böylece bir genom veya ekzomun tamamı saatler yerine dakikalar içinde tamamlanır. Bir varyant, SpliceAI'nin ilgili gen için tanımladığı aralıkta yer alıyorsa (skorun hangi gene ait olduğunu belirleyen kural) o genin içinde kabul edilir ve yalnızca bu genlere ait satırlar tutulur.\n\nÇalıştır, ikisinden biri seçilince açılır: Tüm genler her varyantı skorlar; Gen listesi seç..., yapıştırılan bir listeyi (örneğin bir HPO fenotip listesi veya gen paneli) kabul eder ve listeyi daha sonra kullanmak üzere bir adla kaydedebilir. Son liste yeniden önerilir.\n\nBu özellik, sonuç tablosunun üstündeki gen seti filtresinden farklıdır; o filtre yalnızca tüm varyantların skorlandığı bir çalışmadaki satırları gizler.",
    '"Check compatibility" performs a quick structural check of the VCF text currently loaded or pasted, without running the scoring, so that problems can be found before a potentially long run.\n\nIt verifies that:\n  • a #CHROM header line with the required columns (#CHROM, POS, ID, REF, ALT) is present\n  • at least one variant record follows the header\n  • a sample of records can be parsed: POS is a positive integer, REF contains valid bases and ALT is not empty\n\nIt also reports the number of variant records and whether the chromosome names carry the "chr" prefix.\n\nTypical variant counts (these vary considerably between laboratories, capture kits and filtering):\n  • WES (whole exome): tens of thousands\n  • WGS (whole genome): several million\n  • CES / clinical exome or gene panels: from hundreds to a few thousand':
        '"Uyumluluğu kontrol et", skorlamayı çalıştırmadan, yüklenen veya yapıştırılan VCF metni üzerinde hızlı bir yapısal kontrol yapar; böylece uzun sürebilecek bir çalışmadan önce olası sorunlar saptanabilir.\n\nKontrol edilenler:\n  • gerekli sütunları (#CHROM, POS, ID, REF, ALT) içeren bir #CHROM başlık satırının bulunması\n  • başlıktan sonra en az bir varyant kaydının bulunması\n  • örneklenen kayıtların doğru ayrıştırılabilmesi: POS pozitif bir tam sayıdır, REF geçerli bazlardan oluşur ve ALT boş değildir\n\nAyrıca varyant kaydı sayısı ve kromozom adlarının "chr" önekini taşıyıp taşımadığı bildirilir.\n\nTipik varyant sayıları (laboratuvara, yakalama kitine ve filtrelemeye göre önemli ölçüde değişir):\n  • WES (tüm ekzom): on binlerce\n  • WGS (tüm genom): birkaç milyon\n  • CES / klinik ekzom veya gen panelleri: yüzlerden birkaç bine kadar',
    'There is one field per genome build. A run uses the field for the build selected above, which is detected automatically from the VCF header when the header states it (##reference, or the length of chr1 in its ##contig lines).\n\nThe simplest option is "Download..." next to a field: it downloads the build\'s FASTA from UCSC (about 1 GB compressed, about 3 GB decompressed), decompresses it and creates its .fai index.\n\nAlternatively, download hg19.fa.gz or hg38.fa.gz from hgdownload.soe.ucsc.edu/goldenPath/{build}/bigZips/, decompress it, create a .fai index (for example with `samtools faidx`), and select the .fa file with Browse.\n\nThe reference FASTA is required even when precomputed data are selected, because variants not covered by them are scored live using the reference sequence.':
        'Her genom build\'i için bir alan bulunur. Çalışmada, yukarıda seçilen build\'e ait alan kullanılır; build, VCF başlığında belirtiliyorsa (##reference veya ##contig satırlarındaki chr1 uzunluğu) otomatik olarak belirlenir.\n\nEn kolay yöntem, alanın yanındaki "İndir..." düğmesidir: build\'e ait FASTA dosyasını UCSC\'den indirir (sıkıştırılmış yaklaşık 1 GB, açıldığında yaklaşık 3 GB), açar ve .fai indeksini oluşturur.\n\nAlternatif olarak hg19.fa.gz veya hg38.fa.gz dosyası hgdownload.soe.ucsc.edu/goldenPath/{build}/bigZips/ adresinden indirilip açılabilir, .fai indeksi oluşturulabilir (örneğin `samtools faidx` ile) ve .fa dosyası Gözat ile seçilebilir.\n\nReferans FASTA, precomputed veri seçildiğinde de gereklidir; çünkü bu verilerde bulunmayan varyantlar referans dizisi kullanılarak live skorlanır.',
    'SnpEff adds transcript and coding-position annotation (see "About SnpEff" (?)). SnpEff is a separate Java program, not part of this Python program.\n\nOne-time setup, if needed:\n  "Download SnpEff..." below downloads SnpEff and the RefSeq transcript database for the selected build (GRCh37.p13/GRCh38.p14, which provide NM_ identifiers), with a progress bar as for the reference FASTA.\n  SnpEff requires Java 21 or later. If this computer does not have it, the same download also obtains a private copy of Java (Eclipse Temurin, about 50 MB) in the program\'s "java" folder, used only by this program; any Java already installed is left unchanged.\n\nOnce set up, Run fills in two additional columns.':
        'SnpEff, transkript ve kodlayan dizi konumu anotasyonunu ekler (ayrıntılar "SnpEff Hakkında" (?) bölümündedir). SnpEff, bu programın parçası olmayan ayrı bir Java programıdır.\n\nGerekirse tek seferlik kurulum:\n  Aşağıdaki "SnpEff indir..." düğmesi, SnpEff\'i ve seçilen build için RefSeq transkript veritabanını (NM_ kimlikleri sağlayan GRCh37.p13/GRCh38.p14) referans FASTA\'da olduğu gibi bir ilerleme çubuğuyla indirir.\n  SnpEff, Java 21 veya daha yeni bir sürüm gerektirir. Bu bilgisayarda bulunmuyorsa aynı indirme, yalnızca bu programın kullandığı özel bir Java kopyasını (Eclipse Temurin, yaklaşık 50 MB) programın "java" klasörüne indirir; bilgisayarda kurulu olan Java değiştirilmez.\n\nKurulumdan sonra Çalıştır iki ek sütunu da doldurur.',
    "SnpEff is an open-source variant annotation tool that determines the gene and transcript in which a variant lies and converts its genomic position into coding-sequence notation (e.g. NM_018136.5, c.875A>T).\n\nWhen a variant overlaps more than one transcript, the transcript shown is chosen in this order: the gene's NCBI MANE Select transcript, if available (matched by RefSeq accession, independent of the genome build); otherwise a curated RefSeq transcript (NM_/NR_) in preference to a computationally predicted one (XM_/XR_); otherwise the first transcript listed by SnpEff.\n\nSpliceAI predicts effects on splicing but does not provide transcript or coding notation; SnpEff complements it, so that each row shows both the splicing prediction and the standard transcript and coding notation.\n\nSnpEff is open source under the MIT license and was developed by Pablo Cingolani.":
        "SnpEff, bir varyantın yer aldığı gen ve transkripti belirleyen ve genomik konumunu kodlayan dizi gösterimine (ör. NM_018136.5, c.875A>T) dönüştüren açık kaynaklı bir varyant anotasyon aracıdır.\n\nBir varyant birden fazla transkriptle örtüştüğünde gösterilecek transkript şu sırayla seçilir: varsa genin NCBI MANE Select transkripti (RefSeq erişim numarasıyla eşleştirilir ve genom build'inden bağımsızdır); yoksa hesaplamayla öngörülmüş bir transkript (XM_/XR_) yerine küratörlü bir RefSeq transkripti (NM_/NR_); bunlar da yoksa SnpEff'in ilk sıraladığı transkript.\n\nSpliceAI splicing üzerindeki etkiyi öngörür, ancak transkript veya kodlayan dizi gösterimi sağlamaz; SnpEff bu bilgiyi tamamlar ve böylece her satırda hem splicing tahmini hem de standart transkript ve kodlayan dizi gösterimi yer alır.\n\nSnpEff, MIT lisansıyla sunulan açık kaynaklı bir yazılımdır ve Pablo Cingolani tarafından geliştirilmiştir.",
    "NCBI's MANE Select is a single, expert-curated representative transcript per gene, agreed by NCBI and Ensembl. When a variant overlaps more than one transcript, SnpEff's transcript and coding-position columns report this transcript in preference to the others.\n\nIt is optional: without it, a curated RefSeq transcript (NM_/NR_) is preferred to a predicted one (XM_/XR_), and otherwise the first transcript listed by SnpEff is used.\n\nThe file is small (about 1 MB) and is downloaded from NCBI's FTP site; the latest release is located automatically.":
        "NCBI MANE Select, NCBI ve Ensembl'in ortaklaşa belirlediği, uzmanlarca küratörlüğü yapılmış ve her gen için tek olan temsili transkripttir. Bir varyant birden fazla transkriptle örtüştüğünde SnpEff'in transkript ve kodlayan dizi konumu sütunlarında öncelikle bu transkript bildirilir.\n\nKullanımı isteğe bağlıdır: kullanılmadığında öngörülmüş bir transkript (XM_/XR_) yerine küratörlü bir RefSeq transkripti (NM_/NR_) tercih edilir; bu da yoksa SnpEff'in ilk sıraladığı transkript kullanılır.\n\nDosya küçüktür (yaklaşık 1 MB) ve NCBI'nin FTP sitesinden indirilir; en güncel sürüm otomatik olarak bulunur.",
    'SpliceAI Variant Scoring':
        'SpliceAI Variant Scoring',
    'VCF scoring':
        'VCF skorlama',
    'SpliceAI ready':
        'SpliceAI hazır',
    'SpliceAI is installed. Click to see its version and location, and whether all its files are there.':
        'SpliceAI kurulu. Sürümünü, konumunu ve tüm dosyalarının yerinde olup olmadığını görmek için tıklayın.',
    "Live scoring can't run without it. Click to install SpliceAI.":
        "Live skorlama onsuz çalışamaz. SpliceAI'yi kurmak için tıklayın.",
    'Building reference DB':
        'Reference DB hesaplanıyor',
    ' (paused)':
        ' (duraklatıldı)',
    'Reference DB: problem':
        'Reference DB: sorun var',
    'Reference DB: {builds}':
        'Reference DB: {builds}',
    'Reference DB: checking...':
        'Reference DB: kontrol ediliyor...',
    'No reference DB':
        'Reference DB yok',
    'Click for details, to build it or import one.':
        'Ayrıntılar, hesaplama ya da içe aktarma için tıklayın.',
    '(version unknown)':
        '(sürüm bilinmiyor)',
    'Location: {path}':
        'Konum: {path}',
    'Model files (spliceai1-5.h5): {found} of {n} found':
        'Model dosyaları (spliceai1-5.h5): {n} dosyadan {found} tanesi bulundu',
    'Gene annotation files (grch37, grch38): {found} of {n} found':
        'Gen anotasyon dosyaları (grch37, grch38): {n} dosyadan {found} tanesi bulundu',
    'SpliceAI check':
        'SpliceAI kontrolü',
    "Missing: {files}. Live scoring won't work until SpliceAI is reinstalled.":
        'Eksik: {files}. SpliceAI yeniden kurulana kadar live skorlama çalışmaz.',
    'Open SpliceAI setup...':
        'SpliceAI kurulumunu aç...',
    'Everything live scoring needs is in place.':
        'Live skorlamanın ihtiyaç duyduğu her şey yerinde.',
    '(This program is tested with SpliceAI {version}.)':
        '(Bu program SpliceAI {version} ile test edilmiştir.)',
    'Input VCF':
        'Girdi VCF',
    'Paste VCF content here, use "Load file..." below, or drag a .vcf / .vcf.gz file anywhere onto this window.':
        'VCF içeriğini buraya yapıştırın, aşağıdaki "Dosya yükle..."yi kullanın ya da bir .vcf / .vcf.gz dosyasını bu pencerenin herhangi bir yerine sürükleyin.',
    'Load file...':
        'Dosya yükle...',
    'Open a VCF file -- or drag one from Explorer anywhere onto this window.':
        "Bir VCF dosyası açın -- ya da Gezgin'den bu pencerenin herhangi bir yerine sürükleyin.",
    "About this VCF -- shows the currently loaded file's reference build, date, sample(s), and variant count, parsed from its own header.":
        "Bu VCF hakkında -- yüklü dosyanın kendi başlığından okunan referans build'ini, tarihini, örnek(ler)ini ve varyant sayısını gösterir.",
    'Settings':
        'Ayarlar',
    'Mode:':
        'Mod:',
    'Raw vs Masked scoring':
        'Raw ve masked skorlama',
    'Use precomputed data':
        'Precomputed veri kullan',
    "I don't have precomputed data":
        'Precomputed verim yok',
    'Precomputed data':
        'Precomputed veri',
    'Precomputed data folder:':
        'Precomputed veri klasörü:',
    'Skip variants already in precomputed data':
        'Precomputed veride zaten olan varyantları atla',
    'Genes to score:':
        'Skorlanacak genler:',
    'Choose genes...':
        'Gen seç...',
    'Choose gene list...':
        'Gen listesi seç...',
    'Choose one to enable Run.':
        "Çalıştır'ın açılması için birini seçin.",
    'All genes':
        'Tüm genler',
    'Genes to score':
        'Skorlanacak genler',
    '{build} reference FASTA (used for this run):':
        '{build} referans FASTA (bu çalışmada kullanılan):',
    '{build}.fa -- keep its .fa.fai next to it':
        '{build}.fa -- .fa.fai dosyası yanında olsun',
    'The {build} reference FASTA -- Browse..., or drag the .fa file here.':
        "{build} referans FASTA'sı -- Gözat... ya da .fa dosyasını buraya sürükleyin.",
    'Download...':
        'İndir...',
    'Download {build}.fa from UCSC, decompress it and build its index':
        "{build}.fa'yı UCSC'den indirir, açar ve indeksini oluşturur",
    'Reference FASTA':
        'Referans FASTA',
    'Annotate with SnpEff':
        'SnpEff ile anotasyon',
    'On by default: SnpEff comes with the program. Turning it off drops the coding position and transcript (NM) columns, and makes a run a little faster.':
        'Varsayılan olarak açık: SnpEff programla birlikte gelir. Kapatmak kodlayan konum ve transkript (NM) sütunlarını kaldırır ve çalışmayı biraz hızlandırır.',
    'SnpEff (always on): gives the coding position and transcript (NM) columns.':
        'SnpEff (her zaman açık): kodlayan konum ve transkript (NM) sütunlarını verir.',
    'MANE Select':
        'MANE Select',
    "On by default: where a gene has a MANE Select transcript, that NM_ is the one reported. Off: SnpEff's transcript order decides, curated (NM_) ones first.":
        "Varsayılan olarak açık: bir genin MANE Select transkripti varsa bildirilen NM_ odur. Kapalıysa SnpEff'in transkript sırası belirler, önce küratörlü (NM_) olanlar.",
    'SnpEff setup':
        'SnpEff kurulumu',
    'SnpEff install folder:':
        'SnpEff kurulum klasörü:',
    'Download SnpEff...':
        'SnpEff indir...',
    'MANE Select folder:':
        'MANE Select klasörü:',
    'Optional: improves the choice of transcript.':
        'İsteğe bağlı: transkript seçimini iyileştirir.',
    'Select a MANE summary file you already have (MANE.GRCh38.vX.X.summary.txt.gz)':
        'Elinizdeki bir MANE özet dosyasını seçin (MANE.GRCh38.vX.X.summary.txt.gz)',
    'Download MANE Select...':
        'MANE Select indir...',
    'MANE Select data':
        'MANE Select verisi',
    'SnpEff wasn\'t found, so results have no transcript (NM) or coding position. It comes with the program: reinstalling restores it, or use "Download SnpEff..." below.':
        'SnpEff bulunamadı; bu nedenle sonuçlarda transkript (NM) ve kodlayan konum yok. SnpEff programla birlikte gelir: yeniden kurulum onu geri getirir ya da aşağıdaki "SnpEff indir..." kullanılabilir.',
    'Pause':
        'Duraklat',
    'Pause the calculation after the current variant; click again to resume.':
        'Hesaplamayı o anki varyanttan sonra duraklatır; devam etmek için yeniden tıklayın.',
    'End':
        'Bitir',
    'End the calculation now (its results are discarded).':
        'Hesaplamayı şimdi bitirir (sonuçları atılır).',
    'Summary':
        'Özet',
    'Summary statistics for the most recently completed run (all results, not the filtered or paginated view).':
        'Son tamamlanan çalışmaya ait özet istatistikler (filtrelenmiş veya sayfalanmış görünüm değil, tüm sonuçlar).',
    "CPU use: this program's share of the whole processor, the whole computer's (including SnpEff's Java), and the time since the run started.":
        "CPU kullanımı: bu programın tüm işlemcideki payı, tüm bilgisayarınki (SnpEff'in Java'sı dahil) ve çalışmanın başından beri geçen süre.",
    'How to interpret scores':
        'Skorlar nasıl yorumlanır',
    'Quick gene filter:':
        'Hızlı gen filtresi:',
    'gene symbol, e.g. BRCA1':
        'gen sembolü, ör. BRCA1',
    'Add':
        'Ekle',
    'Export visible rows...':
        'Görünen satırları dışa aktar...',
    'Exports only the rows currently shown, following the active gene and threshold filters and the sort order.':
        'Yalnızca görüntülenen satırları, etkin gen ve eşik filtreleri ile sıralamaya uygun olarak dışa aktarır.',
    'Download full results...':
        'Tüm sonuçları indir...',
    'Exports the complete result set from the most recent run, ignoring any active filters.':
        'Son çalışmanın tüm sonuç setini, etkin filtreleri yok sayarak dışa aktarır.',
    'Open calculated data...':
        'Hesaplanmış veriyi aç...',
    'Loads a previously saved results file (CLI or GUI output) directly into the table, skipping scoring.':
        'Daha önce kaydedilmiş bir sonuç dosyasını (CLI veya GUI çıktısı) skorlama yapmadan doğrudan tabloya yükler.',
    'Clear results':
        'Sonuçları temizle',
    'Import gene list...':
        'Gen listesi içe aktar...',
    'Clear gene-set filter':
        'Gen seti filtresini temizle',
    '0 rows':
        '0 satır',
    'Analyse transcript effect...':
        'Transkript etkisini analiz et...',
    'Sends the selected variant to the Sequence to cDNA tab: the sample is built from the reference with this variant only, and its predicted transcripts are compared with the reference and the WT prediction. Up to {n} selected variants are run one after another.':
        'Seçili varyantı Sequence to cDNA sekmesine gönderir: örnek, yalnızca bu varyantın uygulandığı referanstan oluşturulur ve tahmin edilen transkriptleri referans ve WT tahminiyle karşılaştırılır. Seçili en fazla {n} varyant arka arkaya çalıştırılır.',
    '< Previous':
        '< Önceki',
    'Page 1 of 1':
        'Sayfa 1/1',
    'Next >':
        'Sonraki >',
    'Calculation still running':
        'Hesaplama hâlâ sürüyor',
    'The calculation is not finished yet. Are you sure you want to exit?\n\nThe run will be stopped and its results discarded.':
        'Hesaplama henüz bitmedi. Çıkmak istediğinizden emin misiniz?\n\nÇalışma durdurulacak ve sonuçları atılacak.',
    'Stopping the calculation...':
        'Hesaplama durduruluyor...',
    'Reference database build running':
        'Referans veritabanı hesaplanıyor',
    'A reference database build is running. Close the program anyway?\n\nThe build stops now and continues from the last finished gene the next time you start the program.':
        'Bir referans veritabanı hesaplaması sürüyor. Program yine de kapatılsın mı?\n\nHesaplama şimdi durur ve programı bir sonraki açışınızda son biten genden devam eder.',
    'Load VCF file':
        'VCF dosyası yükle',
    'VCF files (*.vcf *.vcf.gz *.vcf.bgz *.bgz);;All files (*)':
        'VCF dosyaları (*.vcf *.vcf.gz *.vcf.bgz *.bgz);;Tüm dosyalar (*)',
    "{name} can't be loaded here.\n\nDrag in a VCF (.vcf, .vcf.gz) to score it, a reference FASTA (.fa, .fasta, .fna) for its build's row, or a MANE summary file (MANE.GRCh38.vX.X.summary.txt.gz).":
        '{name} buraya yüklenemez.\n\nSkorlamak için bir VCF (.vcf, .vcf.gz), kendi build satırı için bir referans FASTA (.fa, .fasta, .fna) ya da bir MANE özet dosyası (MANE.GRCh38.vX.X.summary.txt.gz) sürükleyin.',
    'Failed to load file':
        'Dosya yüklenemedi',
    'Loaded: {name}':
        'Yüklendi: {name}',
    'build not stated in the VCF -- make sure this is right':
        "build VCF'te belirtilmemiş -- doğru olduğundan emin olun",
    'from the VCF ({evidence})':
        "VCF'ten ({evidence})",
    'the VCF says {build} ({evidence})':
        'VCF {build} diyor ({evidence})',
    'Genome build: {build}':
        "Genom build'i: {build}",
    "Chosen from the VCF's own header. To set it by hand, open Advanced settings.":
        "VCF'in kendi başlığından seçildi. Elle ayarlamak için Gelişmiş ayarlar'ı açın.",
    '{build} reference FASTA:':
        '{build} referans FASTA:',
    'Hide advanced settings':
        'Gelişmiş ayarları gizle',
    'Advanced settings...':
        'Gelişmiş ayarlar...',
    "the {build} reference genome hasn't been found":
        '{build} referans genomu bulunamadı',
    ' -- open Advanced settings.':
        " -- Gelişmiş ayarlar'ı açın.",
    ', found in {folder}':
        ', bulunduğu yer: {folder}',
    'reference genome {builds}':
        'referans genom {builds}',
    ' and ':
        ' ve ',
    'Ready: {what}{where}.':
        'Hazır: {what}{where}.',
    'No VCF is currently loaded.':
        'Şu anda yüklü VCF yok.',
    'Reference (from header): {ref}':
        'Referans (başlıktan): {ref}',
    'The file is {build}, but {selected} is currently selected.':
        'Dosya {build}, ama şu anda {selected} seçili.',
    'Genome build: not stated in the file (no ##reference line or ##contig lengths)':
        "Genom build'i: dosyada belirtilmemiş (##reference satırı ya da ##contig uzunlukları yok)",
    'File date: {date}':
        'Dosya tarihi: {date}',
    'Sample':
        'Örnek',
    'Samples':
        'Örnekler',
    'Variant records: {n:,}':
        'Varyant kaydı: {n:,}',
    'Source: {source}':
        'Kaynak: {source}',
    'Contigs: {n}':
        'Contig sayısı: {n}',
    'About this VCF':
        'Bu VCF hakkında',
    'PASS':
        'GEÇTİ',
    'FAIL':
        'KALDI',
    'Gene-set filter: {n} gene(s) active':
        'Gen seti filtresi: {n} gen etkin',
    'Score only these genes':
        'Yalnızca bu genleri skorla',
    '{genes} and {n} more':
        '{genes} ve {n} tane daha',
    '{n} gene(s): {genes}':
        '{n} gen: {genes}',
    'all genes':
        'tüm genler',
    'Select precomputed-scores folder':
        'Precomputed skor klasörünü seçin',
    'Select SnpEff install folder':
        'SnpEff kurulum klasörünü seçin',
    'Select the {build} reference FASTA':
        "{build} referans FASTA'sını seçin",
    'FASTA files (*.fa *.fasta *.fna);;All files (*)':
        'FASTA dosyaları (*.fa *.fasta *.fna);;Tüm dosyalar (*)',
    'Different build?':
        'Farklı build?',
    '{name} looks like an {named} file, but this is the {build} field.\n\nUse it for {build} anyway?':
        '{name} bir {named} dosyasına benziyor, ama bu {build} alanı.\n\nYine de {build} için kullanılsın mı?',
    'load a VCF':
        'bir VCF yükleyin',
    'set the {build} reference FASTA':
        "{build} referans FASTA'sını ayarlayın",
    'choose a precomputed-data option':
        'bir precomputed veri seçeneği seçin',
    'choose All genes or a gene list':
        'Tüm genler veya bir gen listesi seçin',
    'set the Precomputed data folder':
        'Precomputed veri klasörünü ayarlayın',
    'Before you can run: {what}.':
        'Çalıştırmadan önce: {what}.',
    "Java for SnpEff isn't installed":
        'SnpEff için Java kurulu değil',
    "SnpEff isn't installed":
        'SnpEff kurulu değil',
    'SnpEff database for {build} not installed':
        '{build} için SnpEff veritabanı kurulu değil',
    '✓ Installed':
        '✓ Kurulu',
    'click "Download SnpEff..."':
        '"SnpEff indir..."e tıklayın',
    'Found (MANE v{version})':
        'Bulundu (MANE v{version})',
    'Not found':
        'Bulunamadı',
    'No MANE.GRCh38.vX.X.summary.txt.gz in {folder}':
        '{folder} içinde MANE.GRCh38.vX.X.summary.txt.gz yok',
    'Select the MANE Select summary file':
        'MANE Select özet dosyasını seçin',
    'MANE summary (MANE.GRCh38.v*.summary.txt.gz);;All files (*)':
        'MANE özeti (MANE.GRCh38.v*.summary.txt.gz);;Tüm dosyalar (*)',
    'Not a MANE summary file':
        'MANE özet dosyası değil',
    "{name} isn't a MANE summary file.\n\nExpected a name like MANE.GRCh38.v1.5.summary.txt.gz (from ftp.ncbi.nlm.nih.gov/refseq/MANE).":
        '{name} bir MANE özet dosyası değil.\n\nMANE.GRCh38.v1.5.summary.txt.gz gibi bir ad bekleniyordu (ftp.ncbi.nlm.nih.gov/refseq/MANE adresinden).',
    'Remove {gene} from the filter':
        '{gene} genini filtreden çıkar',
    'Showing 0 of 0 rows':
        '0/0 satır gösteriliyor',
    'Showing {first:,}–{last:,} of {total:,} rows':
        '{total:,} satırın {first:,}–{last:,} arası gösteriliyor',
    'Page {page} of {pages}':
        'Sayfa {page}/{pages}',
    'No variants found for: {genes}':
        'Şunlar için varyant bulunamadı: {genes}',
    'Copy variant ID ({id})':
        'Varyant kimliğini kopyala ({id})',
    'Select a variant in the results first.':
        'Önce sonuçlarda bir varyant seçin.',
    '{n} variants are selected; at most {most} can be queued.':
        '{n} varyant seçili; en fazla {most} tanesi kuyruğa alınabilir.',
    'reference values from the database':
        'referans değerleri veritabanından',
    'the reference computed live first, {bases:,} bases':
        'önce referans live hesaplanır, {bases:,} baz',
    'Clear results?':
        'Sonuçlar temizlensin mi?',
    'This will discard the current results table. Continue?':
        'Mevcut sonuç tablosu atılacak. Devam edilsin mi?',
    'Genome build mismatch':
        "Genom build'i uyuşmuyor",
    'This VCF is {detected} ({evidence}), but {selected} is selected.\n\nScoring it as {selected} gives wrong or empty results. Run anyway?':
        'Bu VCF {detected} ({evidence}), ama {selected} seçili.\n\nOnu {selected} olarak skorlamak yanlış ya da boş sonuç verir. Yine de çalıştırılsın mı?',
    'Done -- {n} rows.':
        'Bitti -- {n} satır.',
    'Pipeline error':
        'Hesaplama hatası',
    'Run ended -- no results.':
        'Çalışma bitirildi -- sonuç yok.',
    'Resume':
        'Devam et',
    'End the calculation?':
        'Hesaplama bitirilsin mi?',
    'End this calculation now? The variants scored so far are discarded.':
        'Bu hesaplama şimdi bitirilsin mi? Şu ana kadar skorlanan varyantlar atılır.',
    'Ending the calculation...':
        'Hesaplama bitiriliyor...',
    'CPU: this program {pct:.0f}%':
        'CPU: bu program %{pct:.0f}',
    'whole PC {pct:.0f}%':
        'tüm bilgisayar %{pct:.0f}',
    'unknown':
        'bilinmiyor',
    'Total variants processed: {n:,}':
        'İşlenen toplam varyant: {n:,}',
    'Result rows: {n:,} (can exceed variant count -- one row per matching gene/source)':
        'Sonuç satırı: {n:,} (varyant sayısını aşabilir -- eşleşen her gen/kaynak için bir satır)',
    'Source breakdown (result rows):':
        'Kaynağa göre dağılım (sonuç satırları):',
    'Precomputed lookup: {n:,} ({pct})':
        'Precomputed veriden: {n:,} ({pct})',
    'Live-scored: {n:,} ({pct})':
        'Live skorlanan: {n:,} ({pct})',
    'No gene overlap: {n:,} ({pct})':
        'Gen örtüşmesi yok: {n:,} ({pct})',
    'Skipped, not scorable (reason in Source column): {n:,} ({pct})':
        'Atlandı, skorlanamıyor (nedeni Source sütununda): {n:,} ({pct})',
    'Score thresholds (max_score >=):':
        'Skor eşikleri (max_score >=):',
    'Settings used for this run:':
        'Bu çalışmada kullanılan ayarlar:',
    'Build: {build}':
        'Build: {build}',
    'Mode: {mode}':
        'Mod: {mode}',
    'Precomputed data: used':
        'Precomputed veri: kullanıldı',
    'Precomputed data: not used':
        'Precomputed veri: kullanılmadı',
    'SnpEff annotation: enabled':
        'SnpEff anotasyonu: açık',
    'SnpEff annotation: disabled':
        'SnpEff anotasyonu: kapalı',
    'Run duration: {time}':
        'Çalışma süresi: {time}',
    'No results yet':
        'Henüz sonuç yok',
    'Run the pipeline first.':
        'Önce hesaplamayı çalıştırın.',
    'Run Summary':
        'Çalışma özeti',
    'Nothing to export':
        'Dışa aktarılacak bir şey yok',
    'Export visible rows':
        'Görünen satırları dışa aktar',
    'CSV files (*.csv);;TSV files (*.tsv)':
        'CSV dosyaları (*.csv);;TSV dosyaları (*.tsv)',
    'Export complete':
        'Dışa aktarma tamamlandı',
    'Wrote {n} rows to {path}':
        '{n} satır {path} dosyasına yazıldı',
    'Nothing to download':
        'İndirilecek bir şey yok',
    'Download full results':
        'Tüm sonuçları indir',
    'Download complete':
        'İndirme tamamlandı',
    'Wrote {n} rows (full, unfiltered) to {path}':
        '{n} satır (tamamı, filtresiz) {path} dosyasına yazıldı',
    'Open calculated data':
        'Hesaplanmış veriyi aç',
    'Results files (*.csv *.tsv);;All files (*)':
        'Sonuç dosyaları (*.csv *.tsv);;Tüm dosyalar (*)',
    'Loading {name}...':
        '{name} yükleniyor...',
    'Could not load file':
        'Dosya yüklenemedi',
    'Loaded {n} rows from {name}.':
        '{name} dosyasından {n} satır yüklendi.',
    # --- mane_download_dialog.py
    'Download MANE Select data':
        'MANE Select verisini indir',
    "Downloads the current NCBI MANE Select summary file (about 1 MB); the latest version is located automatically in NCBI's directory listing. It is used to prefer the canonical transcript when a variant overlaps several transcripts during SnpEff annotation; it is otherwise optional.":
        "Güncel NCBI MANE Select özet dosyasını (yaklaşık 1 MB) indirir; en son sürüm NCBI'nin klasör listesinden otomatik olarak bulunur. SnpEff anotasyonunda bir varyant birden fazla transkriptle örtüştüğünde kanonik transkriptin seçilmesi için kullanılır; bunun dışında isteğe bağlıdır.",
    'Save to:':
        'Kaydedilecek yer:',
    'Downloading...':
        'İndiriliyor...',
    'MANE Select data ready':
        'MANE Select verisi hazır',
    'Downloaded:\n{path}':
        'İndirildi:\n{path}',
    # --- refdb.py
    'The {build} reference score database in {folder}':
        '{folder} içindeki {build} referans skor veritabanı',
    '{name} has format {format}, and this program reads format {ours}. It was made by another version of the program; build it again.':
        '{name} {format} biçiminde, bu program ise {ours} biçimini okuyor. Programın başka bir sürümüyle yapılmış; yeniden hesaplayın.',
    '{name} is for {other}, not {build}.':
        '{name} {build} için değil, {other} için.',
    "{name} was computed with other SpliceAI model files ({files} differ from the installed ones), so its values wouldn't match this program's own scoring.":
        '{name} başka SpliceAI model dosyalarıyla hesaplanmış ({files} kurulu olanlardan farklı); bu yüzden değerleri bu programın kendi skorlamasıyla uyuşmaz.',
    '{name} failed its validation on {date}':
        '{name}, {date} tarihindeki doğrulamadan geçemedi',
    '{name} failed its verification on {date}':
        '{name}, {date} tarihindeki kontrolden geçemedi',
    '{what} ({problem}). Build it again, or import a database that passed.':
        '{what} ({problem}). Yeniden hesaplayın ya da doğrulamadan geçmiş bir veritabanı içe aktarın.',
    '{what}. Build it again, or import a database that passed.':
        '{what}. Yeniden hesaplayın ya da doğrulamadan geçmiş bir veritabanı içe aktarın.',
    '{name} is incomplete: {file} is missing.':
        '{name} eksik: {file} yok.',
    "{name} is damaged: {file} isn't the file its manifest describes (size or SHA-256 differ). Build it again, or import an intact copy.":
        "{name} bozuk: {file}, manifest'te tarif edilen dosya değil (boyutu ya da SHA-256'sı farklı). Yeniden hesaplayın ya da sağlam bir kopyasını içe aktarın.",
    '{folder} holds no finished reference score database (no manifest.json).':
        '{folder} içinde tamamlanmış bir referans skor veritabanı yok (manifest.json yok).',
    '{folder} holds a database for an unknown genome build ({build}).':
        "{folder} içindeki veritabanı bilinmeyen bir genom build'i için ({build}).",
    '{folder} is already the {build} database this program uses.':
        '{folder} zaten bu programın kullandığı {build} veritabanı.',
    '{folder} is incomplete: {files} missing.':
        '{folder} eksik: {files} yok.',
    # --- refdb_wizard.py
    'The Sequence to cDNA tab requires SpliceAI\'s acceptor and donor probabilities at every base of the reference genome\'s genes. As predicting them takes SpliceAI several minutes per gene, they are computed once and stored in a reference score database: one per genome build (hg19, hg38), about 1.1 GB each, in SpliceAI_reference_data\\spliceai-refdb in the user folder, next to the reference genomes.\n\nWhen the program is installed from the USB drive, Setup copies the database. Otherwise, it can be obtained in two ways:\n  • From a folder: a database folder obtained separately (for example spliceai-refdb\\hg19 on the USB drive) is installed after being checked.\n  • Compute it on this computer: SpliceAI computes the database from the reference genome. All genes take many hours per build (about 9 on a 6-core computer); meanwhile the program remains usable (the VCF tab is not affected), and the computation can be paused or the program closed, as it resumes at the next start.\n\nBefore use, each database is checked: it must be for the selected build, must have been produced with the same SpliceAI model files as those installed, and every file must match its checksum. A database computed on this computer is also compared with a small validation set of known values and is not used if it does not match.\n\nGenes not yet in the database are computed live when needed (a minute or two each), cached and marked "computed live".\n\nAs the database is derived from SpliceAI\'s models, SpliceAI\'s terms apply to it: academic and non-commercial use only.':
        'Sequence to cDNA sekmesi, referans genomdaki genlerin her bazı için SpliceAI\'nin akseptör ve donör olasılıklarına ihtiyaç duyar. Bu olasılıkların hesaplanması gen başına birkaç dakika sürdüğünden bir kez hesaplanarak referans skor veritabanında saklanır: her genom build\'i (hg19, hg38) için bir veritabanı, her biri yaklaşık 1,1 GB, kullanıcı klasöründe referans genomların yanında SpliceAI_reference_data\\spliceai-refdb içinde.\n\nProgram USB bellekten kurulduğunda veritabanı Setup tarafından kopyalanır. Aksi halde iki yolla edinilebilir:\n  • Klasörden: ayrıca edinilmiş bir veritabanı klasörü (örneğin USB bellekteki spliceai-refdb\\hg19) kontrol edildikten sonra kurulur.\n  • Bu bilgisayarda hesapla: SpliceAI veritabanını referans genomdan hesaplar. Tüm genler için işlem build başına saatler sürer (6 çekirdekli bir bilgisayarda yaklaşık 9 saat); bu süre boyunca program kullanılabilir (VCF sekmesi etkilenmez), hesaplama duraklatılabilir ve program kapatılabilir; işlem bir sonraki açılışta kaldığı yerden devam eder.\n\nHer veritabanı kullanılmadan önce kontrol edilir: seçilen build\'e ait olmalı, kurulu olanlarla aynı SpliceAI model dosyalarıyla üretilmiş olmalı ve her dosya kendi sağlama toplamıyla eşleşmelidir. Bu bilgisayarda hesaplanan bir veritabanı ayrıca bilinen değerlerden oluşan küçük bir doğrulama setiyle karşılaştırılır; eşleşmezse kullanılmaz.\n\nHenüz veritabanında bulunmayan genler gerektiğinde live hesaplanır (her biri bir iki dakika), önbellekte saklanır ve "computed live" olarak işaretlenir.\n\nVeritabanı SpliceAI modellerinden türetildiğinden SpliceAI\'nin kullanım koşulları bu veritabanı için de geçerlidir: yalnızca akademik ve ticari olmayan kullanım.',
    'CPU threads: the number of processor threads the computation may use. "All" is fastest; fewer threads leave more of the computer free for other work (the computation already runs at low priority, so the program and other applications remain responsive in either case).\n\nKeep this computer awake: prevents Windows from going to sleep while the computation runs (the screen may still turn off). This is released when the computation finishes, is paused or stopped, or the program closes.\n\nA Windows Update restart interrupts the computation. No data are lost: it resumes from the last completed gene the next time the program starts.':
        'CPU iş parçacığı: hesaplamanın kullanabileceği işlemci iş parçacığı sayısıdır. "Tümü" en hızlı seçenektir; daha az iş parçacığı, bilgisayarın daha büyük bir kısmını diğer işlere bırakır (hesaplama zaten düşük öncelikle çalıştığından program ve diğer uygulamalar her durumda yanıt vermeye devam eder).\n\nBu bilgisayarı uyanık tut: hesaplama sürerken Windows\'un uyku moduna geçmesini engeller (ekran yine de kapanabilir). Bu ayar hesaplama bittiğinde, duraklatıldığında veya durdurulduğunda ya da program kapatıldığında kaldırılır.\n\nWindows Update kaynaklı bir yeniden başlatma hesaplamayı kesintiye uğratır. Bu durumda veri kaybı olmaz; hesaplama, program bir sonraki açıldığında son tamamlanan genden devam eder.',
    'Genes not (yet) in the database are computed live when an analysis needs them, using the same method as the database, and stored here so that each is computed only once. They are marked "computed live". Clearing the cache only means that they will be computed again.':
        'Veritabanında (henüz) bulunmayan genler, bir analiz için gerektiğinde veritabanıyla aynı yöntemle live hesaplanır ve her biri yalnızca bir kez hesaplanacak şekilde burada saklanır. Bu genler "computed live" olarak işaretlenir. Önbelleğin temizlenmesi yalnızca bu genlerin yeniden hesaplanmasını gerektirir.',
    "The validation set contains SpliceAI's values for about 50 genes per build, computed when this program's own databases were built and checked. After a computation on this computer, the genes it shares with the validation set are compared with these values; a database that does not match is not used. The validation set is a separate file of a few MB (it may not be included with the program for licensing reasons); it can be selected here if available.":
        'Doğrulama seti, build başına yaklaşık 50 gen için SpliceAI değerlerini içerir; bu değerler programın kendi veritabanları üretilip kontrol edilirken hesaplanmıştır. Bu bilgisayarda bir hesaplama tamamlandığında, doğrulama setinde de bulunan genler bu değerlerle karşılaştırılır; eşleşmeyen bir veritabanı kullanılmaz. Doğrulama seti birkaç MB büyüklüğünde ayrı bir dosyadır (lisans nedeniyle programla birlikte sunulmayabilir); mevcutsa buradan seçilebilir.',
    "The {build} database in {folder} can't be read: {error}":
        '{folder} içindeki {build} veritabanı okunamıyor: {error}',
    'checking...':
        'kontrol ediliyor...',
    'not on this computer yet':
        'henüz bu bilgisayarda yok',
    'not used -- {why}':
        'kullanılmıyor -- {why}',
    '{n:,} genes (all)':
        '{n:,} gen (tümü)',
    '{n:,} genes':
        '{n:,} gen',
    'updated {date}':
        'güncelleme {date}',
    'verified':
        'doğrulandı',
    'not verified':
        'doğrulanmadı',
    'ready':
        'hazır',
    "Starting (loading SpliceAI's models)...":
        'Başlıyor (SpliceAI modelleri yükleniyor)...',
    '{done:,} of {total:,} genes':
        '{total:,} genden {done:,} tanesi',
    'overall {pct:.1f}%':
        'toplam %{pct:.1f}',
    'about {time} left':
        'yaklaşık {time} kaldı',
    'Paused. ':
        'Duraklatıldı. ',
    'Stopped. ':
        'Durduruldu. ',
    'No {build} reference FASTA is set (Settings > Advanced settings).':
        "{build} referans FASTA'sı ayarlanmamış (Ayarlar > Gelişmiş ayarlar).",
    'The {build} reference FASTA has no index ({name}.fai) next to it.':
        "{build} referans FASTA'sının yanında indeksi ({name}.fai) yok.",
    "SpliceAI isn't installed (the button at the top of the window).":
        'SpliceAI kurulu değil (pencerenin üstündeki düğme).',
    "The curated RefSeq table for {build} isn't there ({folder}).":
        '{build} için küratörlü RefSeq tablosu yok ({folder}).',
    "The MANE Select summary file isn't there (Settings > Advanced settings).":
        'MANE Select özet dosyası yok (Ayarlar > Gelişmiş ayarlar).',
    "The {build} build can't start: {problem}":
        '{build} hesaplaması başlayamıyor: {problem}',
    'Build stopped':
        'Hesaplama durduruldu',
    'The {build} build was stopped. The genes it finished are kept and used; build again to continue it.':
        '{build} hesaplaması durduruldu. Biten genler saklanır ve kullanılır; devam etmek için yeniden hesaplayın.',
    'The {build} build stopped unexpectedly (exit code {code}). Its finished genes are kept; it continues at the next start. Details: {log}':
        '{build} hesaplaması beklenmedik şekilde durdu (çıkış kodu {code}). Biten genler saklanır; bir sonraki açılışta devam eder. Ayrıntılar: {log}',
    'Build stopped unexpectedly':
        'Hesaplama beklenmedik şekilde durdu',
    'Build finished':
        'Hesaplama bitti',
    'The {build} database is built. It could not be verified: no validation set was found (Reference database > Validation set).':
        '{build} veritabanı hesaplandı. Doğrulanamadı: doğrulama seti bulunamadı (Reference database > Validation set).',
    'The {build} database is built, but verifying it failed: {error}':
        '{build} veritabanı hesaplandı ama doğrulaması başarısız oldu: {error}',
    "The {build} database is built. None of the validation set's genes are in it, so it could not be verified.":
        '{build} veritabanı hesaplandı. Doğrulama setinin hiçbir geni içinde olmadığı için doğrulanamadı.',
    "The {build} database's verification couldn't be recorded: {error}":
        '{build} veritabanının doğrulaması kaydedilemedi: {error}',
    'The {build} database is built and verified: {n} genes match the validation set (largest difference {diff:.6f}).':
        '{build} veritabanı hesaplandı ve doğrulandı: {n} gen doğrulama setiyle eşleşiyor (en büyük fark {diff:.6f}).',
    "The {build} database is built but failed its verification, so it won't be used: ":
        '{build} veritabanı hesaplandı ama doğrulamadan geçemedi, bu yüzden kullanılmayacak: ',
    'From a folder (USB drive)...':
        'Klasörden (USB bellek)...',
    'Install a database folder obtained separately (e.g. spliceai-refdb\\hg19 on the USB drive) -- or drag the folder here.':
        'Ayrıca edinilmiş bir veritabanı klasörünü kurar (örneğin USB bellekteki spliceai-refdb\\hg19) -- ya da klasörü buraya sürükleyin.',
    'Compute it on this computer...':
        'Bu bilgisayarda hesapla...',
    'Compute it here with SpliceAI: many hours per genome build.':
        "SpliceAI ile burada hesaplar: genom build'i başına saatler sürer.",
    'Choose the database folder (the one with manifest.json)':
        "Veritabanı klasörünü seçin (manifest.json'ın bulunduğu)",
    'A database is being imported; wait until it is finished.':
        'Bir veritabanı içe aktarılıyor; bitmesini bekleyin.',
    'Import database':
        'Veritabanını içe aktar',
    'The {build} database is being built right now. Stop the build first, then import.':
        '{build} veritabanı şu anda hesaplanıyor. Önce hesaplamayı durdurun, sonra içe aktarın.',
    'Replace the database?':
        'Veritabanı değiştirilsin mi?',
    'There already is an {build} database. Replace it with the imported one? It is replaced only after the copy is complete and checked.':
        'Zaten bir {build} veritabanı var. İçe aktarılanla değiştirilsin mi? Değişiklik ancak kopyalama tamamlanıp kontrol edildikten sonra yapılır.',
    'Copying and checking {folder}...':
        '{folder} kopyalanıyor ve kontrol ediliyor...',
    'Import failed':
        'İçe aktarma başarısız',
    'Database imported':
        'Veritabanı içe aktarıldı',
    'The database was checked and installed in {folder}.':
        'Veritabanı kontrol edildi ve {folder} klasörüne kuruldu.',
    'No reference score database yet':
        'Henüz referans skor veritabanı yok',
    "This tab compares a sample's sequence with the reference genome using SpliceAI's predictions at every base of a gene. These predictions for the reference are stored in a database on this computer (about 1.1 GB per genome build). The database is provided on the program's USB drive (folder spliceai-refdb; Setup run from the drive copies it) or can be computed on this computer. The VCF tab does not require it.":
        "Bu sekme, bir örneğin dizisini, genin her bazı için SpliceAI tahminlerini kullanarak referans genomla karşılaştırır. Referansa ait bu tahminler bu bilgisayardaki bir veritabanında saklanır (genom build'i başına yaklaşık 1,1 GB). Veritabanı programın USB belleğiyle birlikte sağlanır (spliceai-refdb klasörü; Setup bellekten çalıştırıldığında kopyalanır) veya bu bilgisayarda hesaplanabilir. VCF sekmesi için veritabanı gerekmez.",
    'A Windows Update restart interrupts the build; it then continues at the next start.':
        'Windows Update yeniden başlatması hesaplamayı keser; bir sonraki açılışta kaldığı yerden devam eder.',
    'Stop':
        'Durdur',
    'Stop this build. The genes it finished are kept and used.':
        'Bu hesaplamayı durdurur. Biten genler saklanır ve kullanılır.',
    'Choose...':
        'Seç...',
    'Choose the {build} validation set file (validation-set-{build}.h5)':
        '{build} doğrulama seti dosyasını seçin (validation-set-{build}.h5)',
    'Validation set':
        'Doğrulama seti',
    'Clear cache':
        'Önbelleği temizle',
    'Genes computed live':
        'Live hesaplanan genler',
    '{build} validation set: {path}':
        '{build} doğrulama seti: {path}',
    'none found':
        'bulunamadı',
    'Building {build} -- {n:,} chosen gene(s)':
        '{build} hesaplanıyor -- seçilen {n:,} gen',
    'Building {build} -- all genes':
        '{build} hesaplanıyor -- tüm genler',
    'Genes computed live (cache): {mb:,.1f} MB':
        'Live hesaplanan genler (önbellek): {mb:,.1f} MB',
    'Stop the build?':
        'Hesaplama durdurulsun mu?',
    'Stop this build? The genes it finished are kept and used, and building again continues from there.':
        'Bu hesaplama durdurulsun mu? Biten genler saklanır ve kullanılır; yeniden başlatınca oradan devam eder.',
    'Choose the {build} validation set':
        '{build} doğrulama setini seçin',
    'Validation set (validation-set-{build}.h5);;HDF5 files (*.h5)':
        'Doğrulama seti (validation-set-{build}.h5);;HDF5 dosyaları (*.h5)',
    'Compute the database on this computer':
        'Veritabanını bu bilgisayarda hesapla',
    'Genome build(s):':
        "Genom build'i:",
    'All genes (recommended)':
        'Tüm genler (önerilen)',
    'Only these genes:':
        'Yalnızca bu genler:',
    'none chosen':
        'seçilmedi',
    'A partial build also includes the about 50 genes of the validation set, so that it can be verified. Genes left out are computed live when they are needed.':
        'Kısmi bir hesaplama, doğrulanabilmesi için doğrulama setinin yaklaşık 50 genini de içerir. Dışarıda kalan genler gerektiğinde live hesaplanır.',
    'CPU threads:':
        'CPU iş parçacığı:',
    'All':
        'Tümü',
    'Build options':
        'Hesaplama seçenekleri',
    'Keep this computer awake during the build':
        'Hesaplama sürerken bu bilgisayarı uyanık tut',
    'It runs in the background: the program remains usable, and closing it (or a Windows Update restart) only pauses the computation until the next start. Genes can be used as soon as they are completed. Progress is shown under Reference DB at the top of the window.':
        'Hesaplama arka planda çalışır: program kullanılabilir durumda kalır; programın kapatılması (veya Windows Update kaynaklı bir yeniden başlatma) hesaplamayı yalnızca bir sonraki açılışa kadar duraklatır. Tamamlanan genler hemen kullanılabilir. İlerleme, pencerenin üstündeki Reference DB düğmesinden izlenebilir.',
    'Genes to build':
        'Hesaplanacak genler',
    'Use these genes':
        'Bu genleri kullan',
    '{n:,} gene(s):':
        '{n:,} gen:',
    'Checking this computer':
        'Bu bilgisayar kontrol ediliyor',
    "Reading the gene tables and measuring this computer's speed (about a minute)...":
        'Gen tabloları okunuyor ve bu bilgisayarın hızı ölçülüyor (yaklaşık bir dakika)...',
    'The check failed: {error}':
        'Kontrol başarısız: {error}',
    '{n:,} genes, {bases:,} bases -> about {size}':
        '{n:,} gen, {bases:,} baz -> yaklaşık {size}',
    ', about {time}':
        ', yaklaşık {time}',
    ', time unknown (the speed measurement failed)':
        ', süre bilinmiyor (hız ölçümü başarısız oldu)',
    'Not enough free disk space: {need} needed while building, {free} free on that drive.':
        'Yeterli boş disk alanı yok: hesaplama sırasında {need} gerekiyor, o sürücüde {free} boş.',
    'Free disk space: {free} ({need} needed while building).':
        'Boş disk alanı: {free} (hesaplama sırasında {need} gerekir).',
    'Measured speed: {speed:,} bases/s (a whole build runs somewhat slower; the times above allow for that).':
        'Ölçülen hız: {speed:,} baz/sn (tüm hesaplama biraz daha yavaş çalışır; yukarıdaki süreler bunu hesaba katar).',
    'Reference FASTA, SpliceAI and the gene tables are in place.':
        'Referans FASTA, SpliceAI ve gen tabloları yerinde.',
    'Build the reference score database':
        'Referans skor veritabanını hesapla',
    'Start building':
        'Hesaplamayı başlat',
    # --- reference_note.py
    "<p><b>SpliceAI Variant Scoring is installed and ready.</b></p><p>Scoring variants requires a <b>reference genome FASTA</b> (hg19.fa or hg38.fa, about 3 GB each). Because of their size, these files are not included in Setup; this does not indicate a missing or faulty installation. The program looks for them each time it starts, and this note appears only when they cannot be found.</p><ul><li><b>If the USB drive used for Setup is available</b>, the quickest solution is to run Setup again from it and leave <b>Copy the reference genomes to this PC</b> selected. The program then finds the files by itself.</li><li><b>If the files are already on this computer</b>, place the <b>.fa</b> file and its <b>.fa.fai</b> file together in a folder named <b>SpliceAI_reference_data</b> in the user folder; they will be found at the next start. They can also be selected manually: <b>Advanced settings...</b> at the bottom of <b>Settings</b>, then <b>Browse...</b> on the hg19 or hg38 row. A .fa file can also be dragged onto this window.</li><li><b>Otherwise</b>, <b>Download...</b> on a row downloads that build from UCSC (about 1 GB). On hospital or institutional networks, this download may be refused with a certificate message; copying the files from the USB drive avoids this.</li><li>Illumina's precomputed score files are optional. If they are not available, select <b>I don't have precomputed data</b>; every variant is then scored directly by SpliceAI.</li></ul>":
        "<p><b>SpliceAI Variant Scoring kuruldu ve kullanıma hazır.</b></p><p>Varyantların skorlanabilmesi için bir <b>referans genom FASTA</b> dosyası gereklidir (hg19.fa veya hg38.fa, her biri yaklaşık 3 GB). Bu dosyalar boyutları nedeniyle Setup'a dahil edilmemiştir; bu durum eksik veya hatalı bir kurulum anlamına gelmez. Program her açılışta bu dosyaları arar ve bu not yalnızca dosyalar bulunamadığında gösterilir.</p><ul><li><b>Setup'ın çalıştırıldığı USB bellek mevcutsa</b>, en hızlı çözüm Setup'ı bellekten yeniden çalıştırmak ve <b>Referans genomları bu bilgisayara kopyala</b> seçeneğini işaretli bırakmaktır. Bundan sonra program dosyaları kendiliğinden bulur.</li><li><b>Dosyalar bu bilgisayarda mevcutsa</b>, <b>.fa</b> dosyası ve <b>.fa.fai</b> dosyası kullanıcı klasöründeki <b>SpliceAI_reference_data</b> adlı klasöre birlikte konulduğunda bir sonraki açılışta bulunur. Dosyalar elle de seçilebilir: <b>Ayarlar</b> bölümünün altındaki <b>Gelişmiş ayarlar...</b>, ardından hg19 veya hg38 satırındaki <b>Gözat...</b>. Bir .fa dosyası doğrudan bu pencereye de sürüklenebilir.</li><li><b>Bunlar mümkün değilse</b>, ilgili satırdaki <b>İndir...</b> düğmesi o build'i UCSC'den indirir (yaklaşık 1 GB). Hastane veya kurum ağlarında bu indirme bir sertifika uyarısıyla engellenebilir; dosyaların USB bellekten kopyalanması bu sorunu önler.</li><li>Illumina'nın precomputed skor dosyaları isteğe bağlıdır. Bu dosyalar mevcut değilse <b>Precomputed verim yok</b> seçeneği seçilmelidir; bu durumda her varyant doğrudan SpliceAI ile skorlanır.</li></ul>",
    'Reference data':
        'Referans verisi',
    "Don't show this again":
        'Bunu bir daha gösterme',
    # --- sequence_view.py
    'Sample candidate:':
        'Örnek adayı:',
    'Sample {rank} (relative score {score:.3f}; {evidence})':
        'Sample {rank} (göreli skor {score:.3f}; {evidence})',
    'Copy cDNAs':
        "cDNA'ları kopyala",
    'Copy proteins':
        'Proteinleri kopyala',
    'No predicted transcripts to compare.':
        'Karşılaştırılacak tahmin edilmiş transkript yok.',
    "colours: the reference exon each base's exon matches (grey: novel), underlined: the first base of an exon, red: a gap.":
        'renkler: her bazın ekzonunun eşleştiği referans ekzon (gri: yeni), altı çizili: bir ekzonun ilk bazı, kırmızı: boşluk.',
    'the ATG is lost in this candidate.':
        'bu adayda ATG kaybolmuş.',
    'Protein':
        'Protein',
    'no reference protein.':
        'referans protein yok.',
    'reference {ref} aa, sample {rank}: {sample} aa ({hgvs}); bold red: the first changed residue, orange: frameshifted, black: the new stop.':
        'referans {ref} aa, örnek {rank}: {sample} aa ({hgvs}); kalın kırmızı: değişen ilk kalıntı, turuncu: çerçevesi kaymış, siyah: yeni stop.',
    # --- snpeff_download_dialog.py
    'Downloading Java (SnpEff needs it)':
        'Java indiriliyor (SnpEff için gerekli)',
    'Unpacking Java':
        'Java açılıyor',
    'Downloading SnpEff':
        'SnpEff indiriliyor',
    'Unzipping':
        'Zip açılıyor',
    'Downloading database':
        'Veritabanı indiriliyor',
    'Download SnpEff':
        'SnpEff indir',
    'Install to folder:':
        'Kurulacak klasör:',
    'Downloads SnpEff (about 65 MB) and the RefSeq transcript database for {build} (a few hundred MB) and installs both; depending on the connection, this may take a few minutes. SnpEff requires Java 21 or later; if this computer does not have it, a private copy of Java (about 50 MB) is downloaded automatically.':
        "SnpEff'i (yaklaşık 65 MB) ve {build} için RefSeq transkript veritabanını (birkaç yüz MB) indirip kurar; bağlantı hızına bağlı olarak birkaç dakika sürebilir. SnpEff, Java 21 veya daha yeni bir sürüm gerektirir; bu bilgisayarda bulunmuyorsa programa özel bir Java kopyası (yaklaşık 50 MB) otomatik olarak indirilir.",
    'SnpEff ready':
        'SnpEff hazır',
    'SnpEff and the {build} database are installed:\n{path}':
        'SnpEff ve {build} veritabanı kuruldu:\n{path}',
    # --- specificity.py
    '{label}: {note}':
        '{label}: {note}',
    '{label}: observed: top-1 {top}, max delta {delta:.2f}; random: {changed}/{n} changed, median max delta {median:.2f} (range {low:.2f}-{high:.2f})':
        '{label}: gözlenen: en üstteki transkript {top}, en büyük delta {delta:.2f}; rastgele: {n} denemenin {changed} tanesinde değişti, en büyük deltaların ortancası {median:.2f} (aralık {low:.2f}-{high:.2f})',
    'changed':
        'değişti',
    'unchanged':
        'değişmedi',
    'The {kind} at {pos:,} ({size:,} bp)':
        '{pos:,} konumundaki {kind} ({size:,} bp)',
    'no random replacement for this kind or size of change; left out':
        'bu tür ya da boyuttaki değişiklik için rastgele yedek yok; dışarıda bırakıldı',
    'every random replacement overlapped another change; left out':
        'her rastgele yedek başka bir değişiklikle çakıştı; dışarıda bırakıldı',
    '{label}: the top transcript stays the same for {same} of {n} shifted breakpoints; it changes with {changed} bp':
        '{label}: kaydırılan {n} kırılma noktasının {same} tanesinde en üstteki transkript aynı kalıyor; {changed} bp ile değişiyor',
    '{label}: the top transcript stays the same for {same} of {n} shifted breakpoints':
        '{label}: kaydırılan {n} kırılma noktasının {same} tanesinde en üstteki transkript aynı kalıyor',
    'no breakpoint to shift (a substitution)':
        'kaydırılacak kırılma noktası yok (tek baz değişikliği)',
    # --- spliceai_setup.py
    'installed by Setup, with the program':
        'programla birlikte Setup tarafından kuruldu',
    'installed from within the program':
        'programın içinden kuruldu',
    'installed in this Python environment':
        'bu Python ortamına kurulu',
    # --- spliceai_setup_dialog.py
    '<p><b>SpliceAI is not installed yet.</b></p><p>Setup normally installs it; this occurs when the computer was not connected to the internet during Setup, or when the program is run from source.</p><p>SpliceAI is developed by Illumina and is not included in this program. <b>Download and install</b> downloads the official SpliceAI {version} package (about 16 MB) once from pypi.org and verifies that it is the published file. SpliceAI is licensed for academic and non-commercial use only (source code: PolyForm Strict 1.0.0, models: CC BY-NC 4.0).</p><p><b>Live scoring</b> requires SpliceAI; the rest of the program works without it.</p>':
        "<p><b>SpliceAI henüz kurulu değil.</b></p><p>SpliceAI normalde Setup sırasında kurulur; bu durum, kurulum sırasında internet bağlantısı bulunmadığında veya program kaynak koddan çalıştırıldığında ortaya çıkar.</p><p>SpliceAI, Illumina tarafından geliştirilmiştir ve bu programa dahil değildir. <b>İndir ve kur</b>, resmi SpliceAI {version} paketini (yaklaşık 16 MB) pypi.org'dan bir kez indirir ve yayımlanan dosyayla aynı olduğunu doğrular. SpliceAI yalnızca akademik ve ticari olmayan kullanım için lisanslanmıştır (kaynak kod: PolyForm Strict 1.0.0, modeller: CC BY-NC 4.0).</p><p><b>Live skorlama</b> için SpliceAI gereklidir; programın diğer bölümleri SpliceAI olmadan da çalışır.</p>",
    'SpliceAI Setup Required':
        'SpliceAI kurulumu gerekli',
    'Download and install':
        'İndir ve kur',
    'Install from file...':
        'Dosyadan kur...',
    'Use a copy of {name}, e.g. from a USB drive.':
        '{name} dosyasının bir kopyasını kullanır, örneğin USB bellekteki.',
    'Re-check':
        'Yeniden kontrol et',
    'Continue without SpliceAI':
        'SpliceAI olmadan devam et',
    'Done':
        'Bitti',
    '✓ SpliceAI is installed and importable.':
        '✓ SpliceAI kurulu ve kullanılabilir.',
    '✗ SpliceAI is not currently importable.':
        '✗ SpliceAI şu anda kullanılamıyor.',
    'Installing...':
        'Kuruluyor...',
    'Select the SpliceAI package':
        'SpliceAI paketini seçin',
    'SpliceAI package ({name});;All files (*)':
        'SpliceAI paketi ({name});;Tüm dosyalar (*)',
    "The install reported success, but SpliceAI still isn't importable -- see the log above.":
        'Kurulum başarılı göründü ama SpliceAI hâlâ kullanılamıyor -- yukarıdaki kayda bakın.',
    'Install failed: {message}':
        'Kurulum başarısız: {message}',
    # --- stacked_view.py
    "transcript orientation 5'→3' ({strand}, reverse-complemented)":
        "transkript yönü 5'→3' ({strand}, ters-tamamlayıcı)",
    "genomic orientation ({strand}, 5'→3' points right)":
        "genomik yön ({strand}, 5'→3' sağa bakar)",
    "genomic orientation ({strand}, 5'→3' points left)":
        "genomik yön ({strand}, 5'→3' sola bakar)",
    'Type, size':
        'Tür, boyut',
    'Where':
        'Yer',
    'Nearest sites':
        'En yakın bölgeler',
    'Origin':
        'Kaynak',
    'Sample 1 there':
        'Orada sample 1',
    'origin: {origin}':
        'kaynak: {origin}',
    'duplicate of {what}':
        '{what} kopyası',
    'reference {ref}, sample {sample}':
        'referans {ref}, örnek {sample}',
    'reference: acceptor {acc}, donor {don}':
        'referans: akseptör {acc}, donör {don}',
    'sample: acceptor {acc}, donor {don}':
        'örnek: akseptör {acc}, donör {don}',
    'delta: acceptor {acc}, donor {don}':
        'delta: akseptör {acc}, donör {don}',
    'exon {n}':
        'ekzon {n}',
    'intron':
        'intron',
    'outside':
        'dışında',
    'duplicate of reference exon {n}':
        'referans ekzon {n} kopyası',
    'novel':
        'yeni',
    'Transcript':
        'Transkript',
    'Its exon':
        'Ekzonu',
    'On the haplotype':
        'Haplotip üzerinde',
    'On the reference':
        'Referans üzerinde',
    'Acceptor at start':
        'Başlangıçtaki akseptör',
    'Donor at end':
        'Sondaki donör',
    'vs the reference':
        'Referansa göre',
    'Exon':
        'Ekzon',
    'Event':
        'Olay',
    # --- texts.py
    'deletion':
        'delesyon',
    'insertion':
        'insersiyon',
    'duplication':
        'duplikasyon',
    'inversion':
        'inversiyon',
    'complex':
        'karmaşık değişiklik',
    'snv':
        'SNV',
    'indel':
        'indel',
    'transcript start':
        'transkript başlangıcı',
    'transcript end':
        'transkript bitişi',
    'ATG':
        'ATG',
    'stop codon':
        'stop kodonu',
    'acceptor':
        'akseptör',
    'donor':
        'donör',
    'exon':
        'ekzon',
    # --- transcript_effect_dialog.py
    "Only this variant is applied: the sample is the reference with this one change. The VCF's other variants are not applied.":
        "Yalnızca bu varyant uygulanır: örnek, bu tek değişikliğin yapıldığı referanstır. VCF'teki diğer varyantlar uygulanmaz.",
    'Their phase is unknown: the ones ticked are applied together with this variant on one haplotype, as if all were on the same allele.':
        'Fazları bilinmiyor: işaretlenenler bu varyantla birlikte tek bir haplotipe, hepsi aynı alel üzerindeymiş gibi uygulanır.',
    "The VCF's results have no other scored variant in {gene}.":
        'VCF sonuçlarında {gene} içinde skorlanmış başka varyant yok.',
    'The Sequence-to-cDNA tab will build the sample from the reference ({build}) with the variant given as coordinates, predict its transcripts and compare them with the reference and the WT prediction.':
        'Sequence to cDNA sekmesi, örneği referanstan ({build}) koordinat olarak verilen varyantla oluşturacak, transkriptlerini tahmin edecek ve bunları referans ve WT tahminiyle karşılaştıracak.',
    'The Sequence-to-cDNA tab will build the sample from the reference ({build}) with the {n} variants, one at a time, given as coordinates, predict its transcripts and compare them with the reference and the WT prediction.':
        'Sequence to cDNA sekmesi, örneği referanstan ({build}) koordinat olarak verilen {n} varyantla, birer birer, oluşturacak, transkriptlerini tahmin edecek ve bunları referans ve WT tahminiyle karşılaştıracak.',
    "(can't be written as an event: {error})":
        '(olay olarak yazılamıyor: {error})',
    'Variant:':
        'Varyant:',
    'entered as {line}':
        '{line} olarak girildi',
    'The variant lies in {n} genes: choose the one to analyse.':
        'Varyant {n} genin içinde: analiz edilecek olanı seçin.',
    "Also scored in {genes}, which the Sequence-to-cDNA gene table doesn't have.":
        '{genes} içinde de skorlandı; Sequence to cDNA gen tablosunda bunlar yok.',
    'Canonical transcript:':
        'Kanonik transkript:',
    'Genotype:':
        'Genotip:',
    'Estimated time:':
        'Tahmini süre:',
    'Also apply other variants in this gene':
        'Bu gendeki diğer varyantları da uygula',
    'none: nothing can be predicted':
        'yok: hiçbir şey tahmin edilemez',
    "Each variant is run on its own: only that variant is applied, never the VCF's others. Each result is saved as its own results file in the folder below. At most {n} variants.":
        "Her varyant ayrı çalıştırılır: yalnızca o varyant uygulanır, VCF'teki diğerleri asla uygulanmaz. Her sonuç aşağıdaki klasöre ayrı bir sonuç dosyası olarak kaydedilir. En fazla {n} varyant.",
    'Variant':
        'Varyant',
    'Gene':
        'Gen',
    'Genotype':
        'Genotip',
    'not in the Sequence-to-cDNA gene table':
        'Sequence to cDNA gen tablosunda yok',
    'Save the results in:':
        'Sonuçları şuraya kaydet:',
    'Save the results in':
        'Sonuçları şuraya kaydet',
    'none of its genes ({genes}) is in the Sequence-to-cDNA gene table for {build}.':
        'genlerinden hiçbiri ({genes}) {build} için Sequence to cDNA gen tablosunda yok.',
    '{gene} has no NM_ transcript to predict against.':
        '{gene} geninin tahminde karşılaştırılacak bir NM_ transkripti yok.',
    'it':
        'bu gen',
    "This variant and the {n} ticked below are applied together; the VCF's other variants are not.":
        "Bu varyant ve aşağıda işaretlenen {n} varyant birlikte uygulanır; VCF'teki diğer varyantlar uygulanmaz.",
    'Choose the folder the results files go in.':
        'Sonuç dosyalarının kaydedileceği klasörü seçin.',
    "{folder} isn't a folder.":
        '{folder} bir klasör değil.',
    'about {time} in all':
        'toplam yaklaşık {time}',
    'about {time}':
        'yaklaşık {time}',
    # --- transcripts.py
    "The curated RefSeq transcript table for {build} ({name}) is missing from {folder}. The Sequence to cDNA tab needs it to find genes and their transcripts. Reinstall the program; when running from source, they are in the repository's refseq folder, and installer\\fetch-thirdparty.ps1 fetches them too.":
        '{build} için küratörlü RefSeq transkript tablosu ({name}) {folder} içinde yok. Sequence to cDNA sekmesi genleri ve transkriptlerini bulmak için ona ihtiyaç duyar. Programı yeniden kurun; kaynak koddan çalıştırırken tablolar deponun refseq klasöründedir, installer\\fetch-thirdparty.ps1 de onları indirir.',
    # --- vcf_bridge.py
    'no genotype in the VCF':
        "VCF'te genotip yok",
    '{sample}, the first of {n} samples':
        '{sample}, {n} örneğin ilki',
    '{gt} ({meaning}; {who})':
        '{gt} ({meaning}; {who})',
    '{gt} ({meaning})':
        '{gt} ({meaning})',
    '{text}. {other}':
        '{text}. {other}',
    "the record wasn't found in the VCF text":
        'kayıt VCF metninde bulunamadı',
    'no genotype (the VCF has no GT field)':
        "genotip yok (VCF'te GT alanı yok)",
    'genotype missing':
        'genotip eksik',
    'hemizygous':
        'hemizigot',
    'not carried':
        'taşınmıyor',
    'homozygous':
        'homozigot',
    'Both alleles carry it: no WT transcript is expected from this locus.':
        'İki alel de taşıyor: bu lokustan WT transkript beklenmez.',
    'heterozygous':
        'heterozigot',
    'The other allele is expected to give the WT transcript.':
        'Diğer alelin WT transkripti vermesi beklenir.',
    'heterozygous, another ALT on the other allele':
        'heterozigot, diğer alelde başka bir ALT',
    "The other allele carries another ALT of this record (allele {allele}), so it isn't expected to give the WT transcript.":
        "Diğer alel bu kaydın başka bir ALT'ını taşıyor (alel {allele}); bu yüzden WT transkripti vermesi beklenmez.",
    'not carried by this sample':
        'bu örnek taşımıyor',
    'VCF tab ({mode})':
        'VCF sekmesi ({mode})',
    'VCF tab':
        'VCF sekmesi',
    'cDNA engine':
        'cDNA motoru',
    '{name}: {scores}':
        '{name}: {scores}',
    'the VCF tab has no score for this variant':
        'VCF sekmesinde bu varyant için skor yok',
    'SpliceAI scores no variant that replaces several bases by several others (a delins), and none outside its genes or past its reference; the engine scores it all the same':
        'SpliceAI birkaç bazı birkaç başka bazla değiştiren varyantları (delins) ve genlerinin dışında ya da referansının ötesinde kalanları skorlamaz; motor yine de hepsini skorlar',
    "near the gene's edge: within {n:,} bases of SpliceAI's gene start or end, the VCF tab (get_delta_scores) replaces the sequence beyond the gene with N, while the engine uses the real sequence there (the database's extended mask)":
        "genin ucuna yakın: SpliceAI'nin gen başlangıcı ya da bitişinin {n:,} bazı içinde VCF sekmesi (get_delta_scores) genin ötesindeki diziyi N ile değiştirir, motor ise orada gerçek diziyi kullanır (veritabanının genişletilmiş maskesi)",
    "the VCF tab scored in masked mode: a gain at an annotated splice site and a loss elsewhere count 0 there; the engine's comparison follows SpliceAI's annotation, which can place exon boundaries differently from the canonical transcript shown here":
        "VCF sekmesi masked modda skorladı: anotasyonlu bir splice bölgesindeki kazanç ve başka yerdeki kayıp orada 0 sayılır; motorun karşılaştırması SpliceAI'nin anotasyonunu izler ve bu, ekzon sınırlarını burada gösterilen kanonik transkriptten farklı yerleştirebilir",
    "the VCF tab's value comes from Illumina's precomputed files (their own annotation and rounding), not from a live computation":
        "VCF sekmesinin değeri live bir hesaplamadan değil, Illumina'nın precomputed dosyalarından (kendi anotasyonları ve yuvarlamalarıyla) gelir",
    "other variants were applied on the same haplotype, so the engine's deltas include their effect":
        'aynı haplotipe başka varyantlar da uygulandı; bu yüzden motorun deltaları onların etkisini de içerir',
    'no known reason (the conventions agree here): please report it':
        'bilinen bir neden yok (kurallar burada uyuşuyor): lütfen bildirin',
    # --- vcf_check.py
    'Input is empty -- nothing to check.':
        'Girdi boş -- kontrol edilecek bir şey yok.',
    'No #CHROM header line found.':
        '#CHROM başlık satırı bulunamadı.',
    'Header is missing required column(s): {columns}':
        'Başlıkta gerekli sütun(lar) eksik: {columns}',
    'No variant records found (0 data lines after the header).':
        'Varyant kaydı bulunamadı (başlıktan sonra 0 veri satırı).',
    'line {i}: fewer than 5 tab-separated fields ({n} found)':
        "satır {i}: sekmeyle ayrılmış alan sayısı 5'ten az ({n} bulundu)",
    'line {i}: POS is not a positive integer: {pos}':
        'satır {i}: POS pozitif bir tam sayı değil: {pos}',
    'line {i}: REF has unexpected characters: {ref}':
        'satır {i}: REF beklenmeyen karakterler içeriyor: {ref}',
    'line {i}: ALT is empty':
        'satır {i}: ALT boş',
    'chr-prefixed':
        'chr önekli',
    'bare':
        'öneksiz',
    'unrecognized':
        'tanınmadı',
    ' (sampled first {n} of {total})':
        ' (toplam {total} kaydın ilk {n} tanesi örneklendi)',
    "{bad}/{n} sampled records failed basic parsing{note} -- this doesn't look like well-formed VCF data.":
        'Örneklenen {n} kaydın {bad} tanesi temel ayrıştırmadan geçemedi{note} -- bu, düzgün bir VCF verisine benzemiyor.',
    '{n} variant record(s) found, header OK, chrom naming: {style}.':
        '{n} varyant kaydı bulundu, başlık doğru, kromozom adlandırma: {style}.',
    '{bad}/{n} sampled line(s) had minor issues{note} (see details).':
        'Örneklenen {n} satırın {bad} tanesinde küçük sorunlar var{note} (ayrıntılara bakın).',
    # --- vcf_info.py
    'chr1 length in its ##contig lines':
        '##contig satırlarındaki chr1 uzunluğu',
    # --- view_data.py
    'inserted: {label}':
        'eklenen: {label}',
    '{place} · c.{c}':
        '{place} · c.{c}',
    '{place} · outside the reference transcript':
        '{place} · referans transkriptin dışında',
    'Genomic':
        'Genomik',
    'outside the reference transcript':
        'referans transkriptin dışında',
    'Inserted':
        'Eklenen',
    'Canonical':
        'Kanonik',
    'inverted':
        'inverte',
    'same':
        'aynı',
    'inserted':
        'eklenmiş',
    'deleted':
        'silinmiş',
    '{status} (from the reference)':
        '{status} (referanstan)',
    'Reference':
        'Referans',
    'Status':
        'Durum',
    'Ref. acceptor':
        'Ref. akseptör',
    'Ref. donor':
        'Ref. donör',
    'Sample acceptor':
        'Örnek akseptör',
    'Sample donor':
        'Örnek donör',
    'Delta acceptor':
        'Delta akseptör',
    'Delta donor':
        'Delta donör',
    'Position':
        'Konum',
    'Bases':
        'Bazlar',
    'Scores':
        'Skorlar',
    'Transcripts':
        'Transkriptler',
    'novel exon (its exon {n})':
        'yeni ekzon (kendi ekzonu {n})',
    'duplicate of exon {n}':
        'ekzon {n} kopyası',
    'intron (in reference {region})':
        'intron (referansta {region})',
    # --- viewer_window.py
    'Each column is one base; reference and sample share the columns like a multiple alignment: an insertion in the sample opens columns with a gap in the reference rows, a deletion columns with a gap in the sample\'s, and an inverted stretch is shown in reference order with the sample\'s bases reverse-complemented and marked.\n\nScore rows show SpliceAI\'s probability, at every base, of being an acceptor (blue) or a donor (orange). Zoomed out they are heat strips -- the stronger the colour, the higher the probability, each pixel showing the highest value of the bases it covers, so no strong site is lost --; zoomed in they become bars, with the site threshold as a dashed line; at base level the values are written out. The delta row is the sample\'s value minus the reference\'s (up: higher in the sample).\n\nTranscripts are drawn as exon blocks on a thin intron line whose chevrons point the way the gene is transcribed, coloured by the reference exon each exon matches; grey hatched blocks are exons with no reference counterpart. Beneath each, its own scores at its exonic bases. Sample bases taken from the reference (where the pasted sequence didn\'t reach) are grey; inserted bases have a yellow and inverted ones a lilac background. Hatched score columns have no score (or were scored with N as context near a chromosome end).\n\n"Transcript orientation (5\'→3\')" draws a minus-strand gene reverse-complemented, exon 1 on the left; the line above the view says which orientation is shown. Positions in tooltips are always genomic.\n\nApplied events are bands through every row in their type\'s colour (see the legend), labelled above the ruler with type, size, coordinates and place in the canonical transcript; an insertion is a line in the reference rows and the inserted bases in the sample\'s. "Show origin" outlines, dashed, where inserted bases come from; "Go to origin" goes there. Focus shows only the reference, WT predicted 1 and sample 1 rows. The delta row is scaled to the largest |delta| in view (written at its left); red ticks mark |delta| ≥ 0.1.\n\nMouse wheel: zoom (from the whole gene to single bases); drag: pan; double-click: the whole gene; hover: a short tooltip; click (or Shift + hover): everything about the base or exon, in the panel on the right. Annotated transcripts are drawn with their UTRs thinner than the coding part.':
        'Her sütun bir bazı temsil eder; referans ve örnek, çoklu hizalamada olduğu gibi aynı sütunları paylaşır: örnekteki bir insersiyon referans satırlarında boşluklu sütunlar, bir delesyon ise örnek satırlarında boşluklu sütunlar oluşturur; inversiyon içeren bir bölge referans sırasıyla gösterilir ve örneğin bazları ters-tamamlayıcı olarak işaretlenir.\n\nSkor satırları, SpliceAI\'nin her baz için akseptör (mavi) veya donör (turuncu) olma olasılığını gösterir. Uzaklaştırıldığında bu satırlar ısı şeritleri olarak görünür; renk koyulaştıkça olasılık artar ve her piksel kapsadığı bazların en yüksek değerini gösterdiğinden güçlü bölgeler gözden kaçmaz. Yakınlaştırıldığında çubuklar halinde görünür ve splice bölgesi eşiği kesikli çizgiyle gösterilir; baz düzeyinde değerler sayı olarak yazılır. Delta satırı, örneğin değerinden referansın değerinin çıkarılmasıyla elde edilir (yukarı yön: örnekte daha yüksek).\n\nTranskriptler, ince bir intron çizgisi üzerinde ekzon blokları olarak çizilir; çizgideki oklar genin transkripsiyon yönünü gösterir. Ekzonlar, eşleştikleri referans ekzona göre renklendirilir; gri taralı bloklar referansta karşılığı bulunmayan ekzonlardır. Her transkriptin altında kendi ekzonik bazlarına ait skorlar yer alır. Örnek dizisinde referanstan alınan bazlar (yapıştırılan dizinin kapsamadığı bölgeler) gri, eklenen bazlar sarı, inversiyon bazları ise lila zeminle gösterilir. Taralı skor sütunlarında skor bulunmaz (veya kromozom ucuna yakınlık nedeniyle N bağlamıyla skorlanmıştır).\n\n"Transkript yönü (5\'→3\')", eksi zincirdeki bir geni ters-tamamlayıcı olarak ve 1. ekzon solda olacak şekilde gösterir; görünümün üstündeki satır, gösterilen yönü belirtir. İpuçlarındaki konumlar her zaman genomik konumlardır.\n\nUygulanan olaylar, tüm satırlardan geçen ve olay türünün renginde gösterilen bantlardır (bkz. açıklama); cetvelin üstünde tür, büyüklük, koordinat ve kanonik transkriptteki konumla etiketlenir. Bir insersiyon referans satırlarında bir çizgi, örnek satırlarında ise eklenen bazlar olarak görünür. "Kaynağı göster", eklenen bazların geldiği bölgeyi kesikli çerçeveyle işaretler; "Kaynağa git" bu bölgeye gider. Odak görünümü yalnızca referans, WT predicted 1 ve sample 1 satırlarını gösterir. Delta satırı, görünümdeki en büyük |delta| değerine göre ölçeklenir (bu değer satırın solunda yazılır); kırmızı çentikler |delta| ≥ 0,1 olan konumları işaretler.\n\nFare tekerleği: yakınlaştırır (tüm genden tek baza kadar); sürükleme: kaydırır; çift tıklama: tüm gen; üzerine gelme: kısa bir ipucu; tıklama (veya Shift + üzerine gelme): baz veya ekzonla ilgili her şey sağdaki panelde. Anotasyonlu transkriptlerde UTR\'ler kodlayan bölgeden daha ince çizilir.',
    'Viewer':
        'Görüntüleyici',
    "Transcript orientation (5'→3')":
        "Transkript yönü (5'→3')",
    'Exon schematic':
        'Ekzon şeması',
    'A minus-strand gene drawn reverse-complemented, exon 1 on the left.':
        'Eksi zincirdeki gen ters-tamamlayıcı çizilir, ekzon 1 solda.',
    "A plus-strand gene reads 5'→3' already.":
        "Artı zincirdeki gen zaten 5'→3' okunur.",
    'Stacked view':
        'Katmanlı görünüm',
    'Sequence comparison':
        'Dizi karşılaştırması',
    'Click a base, an exon or an event.':
        'Bir baza, ekzona ya da olaya tıklayın.',
    'Everything about what you clicked shows here (Shift + hover in the stacked view does it too).':
        'Tıkladığınız şeyle ilgili her şey burada görünür (katmanlı görünümde Shift + üzerine gelmek de aynısını yapar).',
    'and {n} more (in the event list below).':
        've {n} tane daha (aşağıdaki olay listesinde).',
    '<b>Applied events</b> (click one to zoom there):':
        '<b>Uygulanan olaylar</b> (oraya yakınlaştırmak için birine tıklayın):',
    'Navigate':
        'Gezin',
    'Show and export':
        'Göster ve dışa aktar',
    'No applied events':
        'Uygulanan olay yok',
    'The applied events; choosing one zooms to it.':
        'Uygulanan olaylar; birini seçmek oraya yakınlaştırır.',
    'Zoom to event':
        'Olaya yakınlaştır',
    'Show origin':
        'Kaynağı göster',
    'Outlines (dashed) where the inserted bases of insertions, duplications, inversions and replacements come from.':
        'İnsersiyon, duplikasyon, inversiyon ve yer değiştirmelerdeki eklenen bazların geldiği yeri kesikli çerçeveyle işaretler.',
    'Go to origin':
        'Kaynağa git',
    "Zooms to where the chosen event's inserted bases come from.":
        'Seçili olayın eklenen bazlarının geldiği yere yakınlaştırır.',
    '◄ Previous':
        '◄ Önceki',
    'The previous difference':
        'Önceki fark',
    'Next ►':
        'Sonraki ►',
    'The next difference':
        'Sonraki fark',
    'Zoom to exon...':
        'Ekzona yakınlaştır...',
    'Whole gene':
        'Tüm gen',
    'Go to: chr7:117,639,961 or c.3718-2477':
        'Git: chr7:117,639,961 ya da c.3718-2477',
    'Focus':
        'Odak',
    "Stacked view: only Reference DNA and scores, the reference transcript, WT predicted 1, the sample's sequence, scores and delta, and sample 1.":
        'Katmanlı görünüm: yalnızca referans DNA ve skorları, referans transkript, WT predicted 1, örneğin dizisi, skorları ve deltası ve sample 1.',
    'Annotated alternatives':
        'Anotasyonlu alternatifler',
    "The gene's other curated transcripts (NM_), from the annotation.":
        'Genin anotasyondaki diğer küratörlü transkriptleri (NM_).',
    'SpliceAI alternatives':
        'SpliceAI alternatifleri',
    'The WT and sample candidates 2 to 5.':
        'WT ve örnek adayları 2-5.',
    'Score tracks in the viewer':
        'Görüntüleyicideki skor izleri',
    '{n:,} columns; showing {first} .. {last}. Reference values: {source}.':
        '{n:,} sütun; gösterilen {first} .. {last}. Referans değerleri: {source}.',
    '{n:,} transcribed position(s) not scored the usual way: ':
        '{n:,} yazılan konum olağan şekilde skorlanmadı: ',
    'Go to':
        'Git',
    "Can't go to {text}: {error}":
        '{text} konumuna gidilemiyor: {error}',
    'Save as PNG':
        'PNG olarak kaydet',
    'PNG (*.png)':
        'PNG (*.png)',
    'Save as SVG':
        'SVG olarak kaydet',
    'SVG (*.svg)':
        'SVG (*.svg)',
    'Export the visible region as CSV':
        'Görünen bölgeyi CSV olarak dışa aktar',
    'CSV (*.csv)':
        'CSV (*.csv)',
    'Open in the viewer':
        'Görüntüleyicide aç',
    # --- worker.py
    'SpliceAI is not installed, so live scoring could not be completed in this run (variants not covered by precomputed data).\n\nTo install it, click the red "SpliceAI not installed" button at the top of the window, then "Download and install" (an internet connection is needed once) or "Install from file...".':
        'SpliceAI kurulu olmadığından bu çalışmada live skorlama tamamlanamadı (precomputed veride bulunmayan varyantlar).\n\nKurulum için pencerenin üstündeki kırmızı "SpliceAI kurulu değil" düğmesine, ardından "İndir ve kur" (bir kez internet bağlantısı gerekir) veya "Dosyadan kur..." seçeneğine tıklayın.',
    'Scores and evidence':
        'Skorlar ve kanıt',
    'Optional checks':
        'İsteğe bağlı kontroller',
    'How to read the result':
        'Sonuç nasıl okunur',
    'Wheel: zoom · drag: pan · double-click: whole gene · Shift + hover: details on the right':
        'Tekerlek: yakınlaştır · sürükle: kaydır · çift tık: tüm gen · Shift + gezdir: ayrıntılar sağda',
    'Compute':
        'Hesapla',
}
