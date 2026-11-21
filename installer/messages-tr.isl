; Setup's own texts in this language (SpliceAI-VariantScoring.iss, [Languages]).
[CustomMessages]
TaskShortcuts=Kısayollar:
TaskDesktopIcon=&Masaüstü kısayolu oluşturulsun
TaskGenomesGroup=Referans genomlar:
TaskGenomes=&Referans genomlar bu bilgisayara kopyalansın (yaklaşık 6 GB)
TaskRefdbGroup=Referans skor veritabanları:
TaskRefdb19=hg&19 referans skor veritabanı kopyalansın (yaklaşık 1,1 GB)
TaskRefdb38=hg&38 referans skor veritabanı kopyalansın (yaklaşık 1,1 GB)
RunNow=%1 başlatılsın
CheckingRefdb=Bu bilgisayardaki %1 referans skor veritabanı denetleniyor...
RefdbKept=- %1: Bu bilgisayarda geçerli bir veritabanı bulunduğundan mevcut veritabanı korunmuştur. Gerektiğinde program içindeki "Klasörden (USB bellek)..." seçeneği ile değiştirilebilir.
RefdbRenamed=- %1: Bu bilgisayardaki veritabanı kullanılamadığından adı %2 olarak değiştirilmiş ve USB bellekteki veritabanı kopyalanmıştır.
RefdbNotRenamed=- %1: Bu bilgisayardaki veritabanı kullanılamadığı ve adı değiştirilemediği için kopyalama yapılmamıştır. Lütfen program içindeki "Klasörden (USB bellek)..." seçeneğini kullanınız.
NoSpaceTitle=Referans verileri için yeterli boş alan bulunmamaktadır
NoSpaceText=Referans genomlarının ve veritabanlarının kopyalanması için yaklaşık %1 GB alan gerekmektedir; bu sürücüde %2 GB boş alan bulunmaktadır.%n%nKurulumu bu dosyaları kopyalamadan sürdürebilir ve dosyaları daha sonra program içinden USB bellekten seçebilirsiniz; bu durumda dosyalar yalnızca USB bellek takılıyken kullanılabilir. Alternatif olarak yer açtıktan sonra kurulumu yeniden başlatabilirsiniz.
NoSpaceInstall=&Kopyalamadan kur
NoSpaceBack=&Geri dön
DownloadTitle=SpliceAI indiriliyor
DownloadText=SpliceAI %1 sürümü pypi.org adresinden indiriliyor (yaklaşık 16 MB).
DownloadFailedTitle=SpliceAI indirilemedi
DownloadFailedText=SpliceAI'nin indirilebilmesi için kurulum sırasında bir kez internet bağlantısı gerekmektedir.%n%nAyrıntılar: %1%n%nLütfen bu bilgisayarın internete bağlı olduğunu kontrol ederek yeniden deneyiniz. İnternet bağlantısı olmadan devam edilmesi durumunda diğer tüm bileşenler kurulur ve program açıldığında SpliceAI'nin indirilmesini önerir.
DownloadRetry=&Yeniden dene
DownloadSkip=&SpliceAI olmadan devam et
UnpackFailed=SpliceAI indirildi ancak açılamadı. Program açıldığında SpliceAI'nin kurulmasını önerecektir.
RefdbNotesTitle=Referans skor veritabanları:

; Inno Setup's own Turkish texts (compiler:Languages\Turkish.isl), in a more formal register.
[Messages]
SelectLanguageTitle=Kurulum dili
SelectLanguageLabel=Lütfen kurulum süresince kullanılacak dili seçiniz.
ClickNext=Devam etmek için Sonraki, kurulumdan çıkmak için İptal düğmesine tıklayınız.
WelcomeLabel2=Bu yardımcı, bilgisayarınıza [name/ver] uygulamasını kuracaktır.%n%nDevam etmeden önce açık olan diğer tüm uygulamaları kapatmanız önerilir.
LicenseLabel=Lütfen devam etmeden önce aşağıdaki önemli bilgileri okuyunuz.
LicenseLabel3=Lütfen aşağıdaki lisans anlaşmasını okuyunuz. Kuruluma devam edebilmek için bu anlaşmayı kabul etmeniz gerekmektedir.
LicenseAccepted=Anlaşmayı kabul &ediyorum
LicenseNotAccepted=Anlaşmayı kabul et&miyorum
WizardSelectTasks=Ek işlemlerin seçimi
SelectTasksDesc=Kurulum sırasında hangi ek işlemler yapılsın?
SelectTasksLabel2=[name] kurulumu sırasında yapılmasını istediğiniz ek işlemleri seçiniz ve Sonraki düğmesine tıklayınız.
WizardReady=Kurulum için hazır
ReadyLabel1=[name] bilgisayarınıza kurulmaya hazırdır.
ReadyLabel2a=Kurulumu başlatmak için Kur düğmesine, ayarları gözden geçirmek veya değiştirmek için Önceki düğmesine tıklayınız.
ReadyLabel2b=Kurulumu başlatmak için Kur düğmesine tıklayınız.
ApplicationsFound=Kurulum tarafından güncellenmesi gereken bazı dosyalar aşağıdaki uygulamalar tarafından kullanılmaktadır. Kurulum yardımcısının bu uygulamaları otomatik olarak kapatmasına izin vermeniz önerilir.
ApplicationsFound2=Kurulum tarafından güncellenmesi gereken bazı dosyalar aşağıdaki uygulamalar tarafından kullanılmaktadır. Kurulum yardımcısının bu uygulamaları otomatik olarak kapatmasına izin vermeniz önerilir. Kurulum tamamlandıktan sonra uygulamalar yeniden başlatılmaya çalışılacaktır.
WizardInstalling=Kurulum sürüyor
InstallingLabel=Lütfen [name] bilgisayarınıza kurulurken bekleyiniz.
FinishedHeadingLabel=[name] kurulumu tamamlandı
FinishedLabelNoIcons=[name] bilgisayarınıza başarıyla kurulmuştur.
FinishedLabel=[name] bilgisayarınıza başarıyla kurulmuştur. Program, oluşturulan kısayollar aracılığıyla başlatılabilir.
ButtonFinish=&Bitir
ClickFinish=Kurulum yardımcısından çıkmak için Bitir düğmesine tıklayınız.
ExitSetupMessage=Kurulum tamamlanmamıştır. Şimdi çıkarsanız program kurulmayacaktır.%n%nKurulumu tamamlamak için kurulum yardımcısını daha sonra yeniden çalıştırabilirsiniz.%n%nKurulum yardımcısından çıkmak istiyor musunuz?
ConfirmUninstall=%1 uygulamasını tüm bileşenleriyle birlikte kaldırmak istediğinizden emin misiniz?
UninstallStatusLabel=Lütfen %1 bilgisayarınızdan kaldırılırken bekleyiniz.
UninstalledAll=%1 bilgisayarınızdan kaldırılmıştır.
UninstalledMost=%1 kaldırılmıştır.%n%nBazı bileşenler kaldırılamamıştır; bunları elle silebilirsiniz.
