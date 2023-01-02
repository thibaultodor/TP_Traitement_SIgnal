#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

#include "Wave.hpp"

double Char2Double(unsigned char val){
    double reel;
    reel = (double)val/127.5;
    return reel-1.0;
}

unsigned char Double2Char(double val){
    val = (val+1.0) * 127.5;
    val = val<0.0?0.0:val;
    val = val>255.0?255.0:val;
    return (unsigned char) floor(val);
}

double * createNote(double freqN,double freqE,int dureeNote){
    int nbEchantillon = dureeNote*(int)freqE;
    auto * data = new double [nbEchantillon];
    for (int i = 0; i < nbEchantillon; ++i) {
        data[i] = sin(2*M_PI*(freqN/freqE)*i);
    }
    return data;
}

void DFT(const double *signal, double *partie_reelle, double *partie_imaginaire, int N){
    double DEUXpi = 2 * M_PI * (1.0/N);
    for (int k = 0; k < N; k++){
        partie_reelle[k] = 0.0;
        partie_imaginaire[k] = 0.0;
        double DeuxPI_K = DEUXpi * k;
        for (int n = 0; n < N; n++){
            partie_reelle[k] += signal[n] * cos(DeuxPI_K * (double)n);
            partie_imaginaire[k] -= signal[n] * sin(DeuxPI_K * (double)n);
        }
    }
}

void IDFT(double *signal, const double *partie_reelle, const double *partie_imaginaire, int N)
{
    double DEUXpi = (2 * M_PI) / N;
    for(int n = 0; n < N; n++)
    {
        double kdpi = DEUXpi * n;
        signal[n] = 0;
        for(int k = 0; k < N; k++){
            signal[n] += (partie_reelle[k] * cos(kdpi * k)) - (partie_imaginaire[k] * sin(kdpi * k));
        }
        signal[n] /= N;
    }
}

unsigned char * gammeChromatique(int secnote,double freqE,std::vector<double> notes){
    int nbEchantillon = (int)freqE*secnote;
    auto * data_char_full = new unsigned char [nbEchantillon*notes.size()];

    for (int i = 0; i < notes.size(); ++i) {
        double * data = createNote(notes[i],freqE,secnote);
        for (int j = 0; j < nbEchantillon; ++j) {
            data_char_full[(i*nbEchantillon)+j] = Double2Char(data[j]);
        }
        delete data;
    }
    return data_char_full;
}

void accord(const std::vector<double>& notes,int secNote,double * reelAccord,double * imgAccord,double freqE){
    int nbEchantillon = secNote * (int)freqE;
    int size = nbEchantillon;

    std::vector<double *> tabReel;
    std::vector<double *> tabImg;
    for (double n:notes) {
        auto * reel = new double [size];
        auto * img = new double [size];
        double * data = createNote(n,freqE,secNote);
        DFT(data,reel,img,size);
        tabReel.emplace_back(reel);
        tabImg.emplace_back(img);
        //delete[] reel;
        //delete[] img;
    }

    for (int i = 0; i < size; ++i) {
        double r = 0;
        double im = 0;
        for (int j = 0; j < notes.size(); ++j) {
            r += tabReel[j][i];
            im += tabImg[j][i];
        }
        reelAccord[i] = r/(int)notes.size();
        imgAccord[i] = im/(int)notes.size();
    }
}

void writeDataFile(const char *s, const double * reel, const double * img, int size){
    FILE* f = fopen(s, "w");
    if( f )
    {
        fprintf(f, "idx,val\n");
        for( int i = 0; i < size; i++ ) {
            double magnitude = sqrt((reel[i] * reel[i]) + (img[i] * img[i]));
            double ft = i;
            fprintf(f, "%f,%f\n", ft, magnitude);
        }
        fclose(f);
    }
    else
        std::cout << "Could not open file !\n";
}

void readDataFile(const char *s,double * reel,double * img){
    std::ifstream input(s);
    std::string delimiter = ",";
    int i = 0;
    for( std::string line; getline( input, line ); )
    {
        if(i!=0){
            size_t pos =0;
            while ((pos = line.find(delimiter)) != std::string::npos) {line.erase(0, pos + delimiter.length());}
            double magnitude = std::stod(line.substr(0, line.find(delimiter)));
            //std::cout<<magnitude<<std::endl;
            reel[i]=magnitude;
            img[i]=magnitude;
        }
        i++;
    }
}

void filtreButterwoth( const double* signal, double* signal_filter, int N , double freq_ech, double freq_c ){
    double alpha = M_PI*freq_c/freq_ech;
    double alpha_2 = alpha*alpha;
    double alpha_3 = alpha_2*alpha;
    int n, k, m;

    double A = 1. + 2.*alpha + 2.*alpha_2 + alpha_3;
    double B = -3. - 2.*alpha + 2.*alpha_2 + 3.*alpha_3;
    double C = 3. - 2.*alpha - 2.*alpha_2 + 3.*alpha_3;
    double D = -1. + 2.*alpha - 2.*alpha_2 + alpha_3;

    double a[4], b[4];
    b[0] = b[3] = alpha_3;
    b[1] = b[2] = 3.*alpha_3;
    a[0] = 0.;
    a[1] = -B/A;
    a[2] = -C/A;
    a[3] = -D/A;

    for( n = 0; n < N; n++ ){
        signal_filter[n] = 0.;
        for( k = 0; k < 4; k++ ){
            m = n-k;
            if( m >= 0 ){
                signal_filter[n] += b[k] * signal[m];
                signal_filter[n] += a[k] * signal_filter[m];
            }
        }
    }
}

void filtreButterwoth( const double* signal, double* signal_filter, int N , double alpha ){
    double alpha_2 = alpha*alpha;
    double alpha_3 = alpha_2*alpha;
    int n, k, m;

    double A = 1. + 2.*alpha + 2.*alpha_2 + alpha_3;
    double B = -3. - 2.*alpha + 2.*alpha_2 + 3.*alpha_3;
    double C = 3. - 2.*alpha - 2.*alpha_2 + 3.*alpha_3;
    double D = -1. + 2.*alpha - 2.*alpha_2 + alpha_3;

    double a[4], b[4];
    b[0] = b[3] = alpha_3;
    b[1] = b[2] = 3.*alpha_3;
    a[0] = 0.;
    a[1] = -B/A;
    a[2] = -C/A;
    a[3] = -D/A;

    for( n = 0; n < N; n++ ){
        signal_filter[n] = 0.;
        for( k = 0; k < 4; k++ ){
            m = n-k;
            if( m >= 0 ){
                signal_filter[n] += b[k] * signal[m];
                signal_filter[n] += a[k] * signal_filter[m];
            }
        }
    }
}

void filtreBessel(const double * signal, double * signal_filter, int size)
{
    double xv[3+1] = {0.0};
    double yv[3+1] = {0.0};

    for (int i=0; i<size; i++)
    {
        xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3];
        xv[3] = signal[i] / double(2.711023309e+02);
        yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3];
        yv[3] = (xv[0] + xv[3]) + 3 * (xv[1] + xv[2])
                + (  0.4226750651 * yv[0]) + ( -1.6550518354 * yv[1])
                + (  2.2028676179 * yv[2]);
        signal_filter[i] = yv[3];
    }
}

Wave::Wave() {            
  //V?rifi que tout est correct pour la taille des types
  checkTypesSize();
  is_data8_allocated  = false; // Tableau non encore allou?
  is_data16_allocated = false; // Tableau non encore allou?
}

Wave::Wave(short* data16,       // Tableau de donn?es
           long int data_nb,    // Nombre de donn?es
           short channels_nb,   // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
           int sampling_freq) { // Fr?quence d'?chantillonnage (en Hertz)(classique en musique : 44100 Hz)

  int i;

  // V?rifi que tout est correct pour la taille des types
  checkTypesSize();

  // Nombre de donn?es
  (*this).data_nb = data_nb;
  // Tableau de don?es lorsque l'on est sur des donn?es 16 bits
  (*this).data16 = new short[data_nb];
  for (i=0; i<data_nb; i++) { 
    (*this).data16[i] = data16[i]; //Recopie en profondeur
  }
  is_data8_allocated  = false; // Tableau non allou?
  is_data16_allocated = true;  // Tableau allou?
  
  InitDescriptor(16, channels_nb, sampling_freq);
}

Wave::Wave(unsigned char* data8, // Tableau de don?es lorsque l'on est sur des donn?es 8 bits
           long int data_nb,     // Nombre de donn?es
           short channels_nb,    // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
           int sampling_freq) {  // Fr?quence d'?chantillonnage (en Hertz)(classique en musique : 44100 Hz)

  int i;

  // V?rifi que tout est correct pour la taille des types
  checkTypesSize();

  // Nombre de donn?es
  (*this).data_nb = data_nb;

  // Tableau de donn?es lorsque l'on est sur des donn?es 16 bits
  (*this).data8 = new unsigned char[data_nb];
  for (i=0; i<data_nb; i++) { 
    (*this).data8[i] = data8[i]; //Recopie en profondeur
  }
  is_data8_allocated  = true;  // Tableau allou?
  is_data16_allocated = false; // Tableau non allou?
  
  InitDescriptor(8, channels_nb, sampling_freq);
}

void Wave::InitDescriptor(int depth,           // Nombre de bits par donn?e (8 ou 16)
                          short channels_nb,   // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
                          int sampling_freq) { // Fr?quence d'?chantillonnage (en Hertz)(classique en musique : 44100 Hz)
           
  // (2 octets) : Nombre de bits par donn?e (8 ou 16)
  (*this).depth = depth;

  // (4 octets) : Constante "RIFF" i.e identification du format
  file_type[0] = 'R'; file_type[1] = 'I'; file_type[2] = 'F'; file_type[3] = 'F';
  // (4 octets) : Identifiant "WAVE"
  file_id[0] = 'W'; file_id[1] = 'A'; file_id[2] = 'V'; file_id[3] = 'E';
  // (4 octets) : Identifiant "fmt "
  chunk_id[0] = 'f'; chunk_id[1] = 'm'; chunk_id[2] = 't'; chunk_id[3] = (unsigned char)32;
  // (4 octets) : Constante "data"
  data_id[0] = 'd'; data_id[1] = 'a'; data_id[2] = 't'; data_id[3] = 'a';          

  // (4 octets) : file_size est le nombre d'octet restant ? lire (i.e = taille du fichier moins 8 octets)
  file_size = 44 + (depth/8) * data_nb - 8;
  // (4 octets) : Nombre d'octets utilis?s pour d?finir en d?tail le chunk
  chunk_size = 16;
  // (2 octets) : Format de fichier (1: PCM,  ...)
  format = 1;
  // (2 octets) : Nombre de canaux (1 pour mono ou 2 pour st?r?o)
  (*this).channels_nb = channels_nb;
  // (4 octets) : Fr?quence d'?chantillonnage (en Hertz)
  (*this).sampling_freq = sampling_freq;
  // (4 octets) : Nombre d'octets par seconde de musique
  (*this).bytes_per_second = sampling_freq*channels_nb*depth/8;
  // (2 octets) : Nombre d'octets par ?chantillon
  (*this).bytes_per_sample = channels_nb*depth/8;
  // (4 octets) : nombre d'octet restant (i.e = taille du fichier moins 44 octets)
  data_size = (depth/8) * data_nb;
}

Wave::~Wave() {            
  if (is_data8_allocated)
    delete[] data8;
  if (is_data16_allocated)
    delete[] data16;
}

void Wave::getData8(unsigned char** data, // Tableau de don?es lorsque l'on est sur des donn?es 8 bits
                    int* size) {          // Taille du tableau
   
  int i;            
  if (!is_data8_allocated) { // Tableau non encore allou?
    cout<<"Wave::getData8: Erreur, les donn?e ne sont pas pr?sente en 8 bits \n";
    exit(-1);
  }

  // Nombre de donn?es
  (*size) = data_nb;
  
  // Allocation du tableau de donn?es 
  (*data) = new unsigned char[data_nb];
  for (i=0; i<data_nb; i++) { 
    (*data)[i] = data8[i]; //Recopie en profondeur
  }
}

void Wave::modifData8(unsigned char* data) { // Tableau de don?es lorsque l'on est sur des donn?es 8 bits
   
  int i;            
  if (!is_data8_allocated) { // Tableau non encore allou?
    cout<<"Wave::setData8: Erreur, le tableau n'est pas allou?\n";
    exit(-1);
  }

  // Remplissage du tableau
  for (i=0; i<data_nb; i++) { 
    data8[i] = data[i];       //Recopie en profondeur
  }
}

void Wave::read(char* fileName) {
  
  int pos = 0;
  int i;
  unsigned char c;
  unsigned char str_tmp[4];
  unsigned char header[44];

  FILE* fd = fopen(fileName, "rb");
  // TEST D'OUVERTURE
  if (fd == NULL) {
    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas ouvrable\n";
    exit(-1);
  }
  // TEST DE LECTURE
  if (fread (header, 1, 44, fd) != 44) {
    cout<<"Wave::read: Erreur, impossible de lire dans le fichier "<<fileName<<" les 44 octets d'entete\n";
    exit(-1);
  }
    
  //for (i=0; i<44; i++) {
  //cout<<(int)header[i]<<endl;
  //}

  // file_type
  for (i=0; i<4; i++, pos++) {
    file_type[i] = header[pos]; // Constante "RIFF" i.e identification du format
    cout<<file_type[i];
  }
  if (!((file_type[0] == 'R') && (file_type[1] == 'I') && (file_type[2] == 'F') && (file_type[3] == 'F'))) {
    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier RIFF\n";
    exit(-1);    
  }
  cout<<endl;

  // file_size
  file_size = 0;
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    file_size +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre d'octets restant ? lire est : "<<file_size<<endl;
  
  // file_id
  for (i=0; i<4; i++, pos++) {
    file_id[i] = header[pos];         // Identifiant "WAVE"
    cout<<file_id[i];
  }
  cout<<endl;
  if (!((file_id[0] == 'W') && (file_id[1] == 'A') && (file_id[2] == 'V') && (file_id[3] == 'E'))) {
    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier WAVE\n";
    exit(-1);    
  }
  
  // chunk_id
  for (i=0; i<4; i++, pos++) {
    chunk_id[i] = header[pos];         // Identifiant "fmt "
    cout<<chunk_id[i];
  }
  cout<<endl;
  if (!((chunk_id[0] == 'f') && (chunk_id[1] == 'm') && (chunk_id[2] == 't') /*&& (chunk_id[3] == ' ')*/)) {
    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier 'fmt'\n";
    exit(-1);    
  }

  // chunk_size
  chunk_size = 0; // Nombre d'octets utilis?s pour d?finir en d?tail le chunk
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    chunk_size +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: chunk_size vaut : "<<chunk_size<<endl;

  // format
  format = 0;     // Format de fichier (1: PCM,  ...)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    format +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le format est (1=PCM) : "<<format<<endl;

  // format
  channels_nb = 0;     // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    channels_nb +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre de canaux est : "<<channels_nb<<endl;

  // sampling_freq
  sampling_freq = 0;     // Fr?quence d'?chantillonnage (en Hertz)
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    sampling_freq +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: La fr?quence d'?chantillonge est : "<<sampling_freq<<endl;

  // bytes_per_second
  bytes_per_second = 0;     // Nombre d'octets par seconde de musique
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    bytes_per_second +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre d'octets par seconde est : "<<bytes_per_second<<endl;

  // bytes_per_sample
  bytes_per_sample = 0;     // Nombre d'octets par ?chantillon
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    bytes_per_sample +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre d'octets par ?chantillon est : "<<bytes_per_sample<<endl;
  if (bytes_per_sample != (bytes_per_second/sampling_freq)) {
    cout<<"Wave::read: bytes_per_sample != (bytes_per_second/sampling_freq)\n";
    exit(-1);
  }

  // depth
  depth = 0;     // Nombre de bits par donn?e (8 ou 16)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    depth +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre de bits par donn?e est : "<<depth<<endl;
  if (depth != (8*bytes_per_sample/channels_nb)) {
    cout<<"Wave::read: depth != (8*bytes_per_sample/channels_nb)\n";
    exit(-1);
  }

  //chunk donn?es
  for (i=0; i<4; i++, pos++) {
    data_id[i] = header[pos]; // Constante "data"
    cout<<data_id[i];
  }
  cout<<endl;
  if (!((data_id[0] == 'd') && (data_id[1] == 'a') && (data_id[2] == 't') && (data_id[3] == 'a'))) {
    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" est corrompu (pas de chaine 'data' dans le header)\n";
    exit(-1);    
  }

  // data_size
  data_size = 0;     // nombre d'octet restant (i.e = taille du fichier moins 44 octets)
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    data_size +=(int)pow(256, i)*(int)header[pos]; //little-endian : increasing numeric significance with increasing memory addresses
  }
  cout<<"Wave::read: Le nombre d'octets restant ? lire est : "<<data_size<<endl;
  if (data_size != (file_size-36)) {
    cout<<"Wave::read: data_size != (file_size-36)\n";
    exit(-1);
  }

  //LECTURE DES DONNEES
  switch (depth) {

    case 8:
      data_nb = data_size;
      data8   = new unsigned char[data_nb]; 
      if (fread (data8, 1, data_nb, fd) != data_nb) { //Les donnees sont sur 8 bits
        cout<<"Wave::read: Erreur, impossible de lire dans le fichier "<<fileName<<" le bon nombre d'octet\n";
        exit(-1);
      }
      is_data8_allocated  = true;  // Tableau non allou?
      is_data16_allocated = false; // Tableau allou?
      break;

    case 16:
      data_nb = data_size/2;
      data16  = new short[data_nb];
      if (fread (data16, 2, data_nb, fd) != data_nb) { //Les donnees sont sur 16 bits
        cout<<"Wave::read: Erreur, impossible de lire dans le fichier "<<fileName<<" le bon nombre d'octet\n";
        exit(-1);
      }
      is_data8_allocated  = false; // Tableau non allou?
      is_data16_allocated = true;  // Tableau allou?
      break;

    default:
      cout<<"Wave::read: Erreur, la profondeur (depth = "<<depth<<") est inconnue\n";
      exit(-1);    
  }

  //FERMETURE DU FICHIER
  fclose(fd);      
}    

void Wave::write(char* fileName) {

  cout <<"Wave::write(char* fileName)\n";
  
  int pos = 0;
  int i;
  unsigned char header[44];

  FILE* fd = fopen(fileName, "wb");
  // TEST D'OUVERTURE
  if (fd == NULL) {
    cout<<"Wave::write: Erreur, le fichier "<<fileName<<" n'est pas ouvrable\n";
    exit(-1);
  }
    
  // file_type : constante "RIFF" i.e identification du format
  for (i=0; i<4; i++, pos++) {
    header[pos] = file_type[i]; // Constante "RIFF" i.e identification du format
  }

  // file_size
  cout<<"file_size = "<< file_size <<endl;
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & file_size)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
//    cout<<"header[pos] = "<<(int)header[pos];
  }
  
  // file_id : identifiant "WAVE"
  for (i=0; i<4; i++, pos++) {
    header[pos] = file_id[i];
  }

  // chunk_id : identifiant "fmt "
  for (i=0; i<4; i++, pos++) {
    header[pos] = chunk_id[i];
  }
  if (!((header[pos-4] == 'f') && (header[pos-3] == 'm') && (header[pos-2] == 't') /*&& (chunk_id[3] == ' ')*/)) {
    cout<<header[pos-4]<<header[pos-3]<<header[pos-2]<<header[pos-1]<<endl;
    cout<<"Wave::write: Erreur, le header est mal remplit pour le chunk_id 'fmt'\n";
    exit(-1);    
  }

  // chunk_size : nombre d'octets utilis?s pour d?finir en d?tail le chunk
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & chunk_size)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // format : format de fichier (1: PCM,  ...)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & format)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // channels_nb : Nombre de canaux (1 pour mono ou 2 pour st?r?o)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & channels_nb)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // sampling_freq : fr?quence d'?chantillonge 
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & sampling_freq)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // bytes_per_second : nombre d'octets par seconde de musique
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & bytes_per_second)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // bytes_per_sample : nombre d'octets par ?chantillon
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & bytes_per_sample)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // depth : nombre de bits par donn?e (8 ou 16)
  for (i=0; i<2; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & depth)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  //chunk donn?es : constante "data"
  for (i=0; i<4; i++, pos++) {
    header[pos] = data_id[i]; // Constante "data"
  }

  // data_size : nombre d'octet restant (i.e = taille du fichier moins 44 octets)
  for (i=0; i<4; i++, pos++) {
    //Bit de poids faible en premier
    header[pos] = (( 255 << (8*i) ) & data_size)>>(8*i); //little-endian : increasing numeric significance with increasing memory addresses
  }

  // ECRITURE DU HEADER
  if (fwrite (header, 1, 44, fd) != 44) { 
    cout<<"Wave::write: Erreur, impossible d'?crire dans le fichier "<<fileName<<" le header\n";
    exit(-1);
  }

  cout<<"Fin ?criture du header\n";
  //LECTURE DES DONNEES
  switch (depth) {
    case 8:
      cout<<"fwrite (data8, \n";
      if (fwrite (data8, 1, data_nb, fd) != data_nb) { //Les donnees sont sur 8 bits
        cout<<"Wave::write: Erreur, impossible d'?crire dans le fichier "<<fileName<<" le bon nombre d'octet\n";
        exit(-1);
      }
      break;
    case 16:
      if (fwrite (data16, 2, data_nb, fd) != data_nb) { //Les donnees sont sur 16 bits
        cout<<"Wave::write: Erreur, impossible d'?crire dans le fichier "<<fileName<<" le bon nombre d'octet\n";
        exit(-1);
      }
      break;
    default:
      cout<<"Wave::write: Erreur, la profondeur (depth = "<<depth<<") est inconnue\n";
      exit(-1);    
  }

  cout<<"Fin ?criture des datas\n";

  //FERMETURE DU FICHIER
  fclose(fd);      
}

//void Wave::read(char* fileName) {
//  
//  int i;
//  unsigned char c;
//  unsigned char str_tmp[4];
//  unsigned char header[44];
//
//  FILE* fd = fopen(fileName, "rb");
//  // TEST D'OUVERTURE
//  if (fd == NULL) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas ouvrable\n";
//    exit(-1);
//  }
//  // TEST DE LECTURE
//  if (fread (header, 1, 44, fd) != 44) {
//    cout<<"Wave::read: Erreur, impossible de lire dans le fichier "<<fileName<<" les 44 octets d'entete\n";
//    exit(-1);
//  }
//  
//  
//  for (i=0; i<44; i++) {
//    cout<<(int)header[i]<<endl;
//  }
//
//
//
//  exit(-1);
//  
//       
//  /*FLUX D'ENTREE POUR LE FICHIER*/	
//  ifstream ifs(fileName,ios::in|ios::binary);
//
//  /*TEST D'OUVERTURE*/
//  if (!ifs.is_open()) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas ouvrable\n";
//    exit(-1);
//  }
//  
//  //LECTURE DU HEADER
//  for (i=0; i<4; i++) {
//    ifs>>file_type[i];       // Constante "RIFF" i.e identification du format
//    cout<<file_type[i];
//  }
//  if (!((file_type[0] == 'R') && (file_type[1] == 'I') && (file_type[2] == 'F') && (file_type[3] == 'F'))) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier RIFF\n";
//    exit(-1);    
//  }
//  cout<<endl;
//  file_size = 0;
//  for (i=0; i<4; i++) {
//    ifs>>str_tmp[i];
//    //Bit de poids faible en premier
//    file_size +=(int)pow(256, i)*(int)str_tmp[i]; //little-endian : increasing numeric significance with increasing memory addresses
//  }
//  cout<<"Wave::read: Le nombre d'octets restant ? lire est : "<<file_size<<endl;
//  
//  for (i=0; i<4; i++) {
//    ifs>>file_id[i];         // Identifiant "WAVE"
//    cout<<file_id[i];
//  }
//  cout<<endl;
//  if (!((file_id[0] == 'W') && (file_id[1] == 'A') && (file_id[2] == 'V') && (file_id[3] == 'E'))) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier WAVE\n";
//    exit(-1);    
//  }
//  for (i=0; i<4; i++) {
//    ifs>>chunk_id[i];         // Identifiant "fmt "
//    cout<<(int)chunk_id[i]<<endl;
//  }
//  cout<<endl;
//  if (!((chunk_id[0] == 'f') && (chunk_id[1] == 'm') && (chunk_id[2] == 't') /*&& (chunk_id[3] == ' ')*/)) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier 'fmt'\n";
//    exit(-1);    
//  }
//
//  // chunk_size
//  chunk_size = 0; // Nombre d'octets utilis?s pour d?finir en d?tail le chunk
//  for (i=0; i<4; i++) {
//    ifs>>str_tmp[i];
//    cout<<(int)str_tmp[i]<<endl;
//    //Bit de poids faible en premier
//    chunk_size +=(int)pow(256, i)*(int)str_tmp[i]; //little-endian : increasing numeric significance with increasing memory addresses
//  }
//    cout<<"Wave::read: chunk_size vaut : "<<chunk_size<<endl;
//
//  // format
//  format = 0;     // Format de fichier (1: PCM,  ...)
//  for (i=0; i<2; i++) {
//    ifs>>str_tmp[i];
//    //Bit de poids faible en premier
//    format +=(int)pow(256, i)*(int)str_tmp[i]; //little-endian : increasing numeric significance with increasing memory addresses
//  }
//  cout<<"Wave::read: Le format est (1=PCM) : "<<format<<endl;
//
//
////  ifs>>channels_nb;          // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
////  cout<<"Wave::read: Le nombre de canaux est : "<<channels_nb<<endl;
////  ifs>>sampling_freq;        // Fr?quence d'?chantillonnage (en Hertz)
////  cout<<"Wave::read: La fr?quence d'?chantillonge est : "<<sampling_freq<<endl;
////  ifs>>bytes_per_second;     // Nombre d'octets par seconde de musique
////  cout<<"Wave::read: Le nombre d'octets par seconde est : "<<bytes_per_second<<endl;
////  ifs>>bytes_per_sample;     // Nombre d'octets par ?chantillon
////  cout<<"Wave::read: Le nombre d'octets par ?chantillon est : "<<bytes_per_sample<<endl;
////  ifs>>depth;                // Nombre de bits par donn?e (8 ou 16)
////  cout<<"Wave::read: Le nombre de bits par donn?e est : "<<depth<<endl;
////*/
////  //chunk donn?es
////  for (i=0; i<4; i++) {
////    ifs>>data_id[i];       // Constante "data"
////    cout<<data_id[i];
////  }
////  if (!((data_id[0] == 'd') && (data_id[1] == 'a') && (data_id[2] == 't') && (data_id[3] == 'a'))) {
////    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" est corrompu (pas de chaine 'data' dans le header)\n";
////    exit(-1);    
////  }
////  ifs>>data_size;       // nombre d'octet restant (i.e = taille du fichier moins 44 octets)
////
////  //LECTURE DES DONNEES
////  switch (depth) {
////    case 8:
////      samples_nb = data_size;
////      data       = new short[samples_nb];
////      char c;
////      for (i=0; i<samples_nb; i++) {
////        ifs>>c; //Lecture de 8 bits
////        data[i] = (short)c;
////      }  
////      break;
////    case 16:
////      samples_nb = data_size/2;
////      data       = new short[samples_nb];
////      for (i=0; i<samples_nb; i++) {
////        ifs>>data[i]; //Lecture de 16 bits
////      }  
////      break;
////    default:
////      cout<<"Wave::read: Erreur, la profondeur (depth = "<<depth<<") est inconnue\n";
////      exit(-1);    
////  }      
//
//  ifs.close();
//}    



//void Wave::read(char* fileName) {
//  
//  int i;
//  char c;
//  char str_tmp[5];
//  str_tmp[4] = '\0';
//       
//  /*FLUX D'ENTREE POUR LE FICHIER*/	
//  FILE* fd = fopen(fileName, "rb");
//
//  /*TEST D'OUVERTURE*/
//  if (fd == NULL) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas ouvrable\n";
//    exit(-1);
//  }
//  //LECTURE DU HEADER
//  for (i=0; i<4; i++) {
//    file_type[i] = (char) fgetc(fd);
//    cout<<file_type[i];
//  }
//  if (!((file_type[0] == 'R') && (file_type[1] == 'I') && (file_type[2] == 'F') && (file_type[3] == 'F'))) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier RIFF\n";
//    exit(-1);    
//  }
//  cout<<endl;
//  for (i=0; i<4; i++) {
//    file_type[i] = (char) fgetc(fd);
//    fscanf(fd, "%c", str_tmp[i]);
//  }
//  file_size = atoi(str_tmp);  // Nombre d'octet restant ? lire (i.e = taille du fichier moins 8 octets)
//  cout<<"Wave::read: nombre d'octet restant ? lire = "<<file_size<<endl;
//  for (i=0; i<4; i++) {
//    fscanf(fd, "%c", file_id[i]);// Identifiant "WAVE"
//    cout<<file_id[i];
//  }
//  cout<<endl;
//  if (!((file_id[0] == 'W') && (file_id[1] == 'A') && (file_id[2] == 'V') && (file_id[3] == 'E'))) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier WAVE\n";
//    exit(-1);    
//  }
//  for (i=0; i<4; i++) {
//    fscanf(fd, "%c", chunk_id[i]);
//    cout<<chunk_id[i];
//  }
//  cout<<endl;
//  if (!((chunk_id[0] == 'f') && (chunk_id[1] == 'm') && (chunk_id[2] == 't') /*&& (chunk_id[3] == ' ')*/)) {
//    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" n'est pas un fichier 'fmt'\n";
//    exit(-1);    
//  }
////  /*TEST*/for (i=0; i<24; i++) {
////  /*TEST*/  ifs>>c;
////  /*TEST*/}
////
/////*  ifs>>chunk_size;           // Nombre d'octets utilis?s pour d?finir en d?tail le chunk
////  ifs>>format;               // Format de fichier (1: PCM,  ...)
////  cout<<"Wave::read: Le format est (1=PCM) : "<<format<<endl;
////  ifs>>channels_nb;          // Nombre de canaux (1 pour mono ou 2 pour st?r?o)
////  cout<<"Wave::read: Le nombre de canaux est : "<<channels_nb<<endl;
////  ifs>>sampling_freq;        // Fr?quence d'?chantillonnage (en Hertz)
////  cout<<"Wave::read: La fr?quence d'?chantillonge est : "<<sampling_freq<<endl;
////  ifs>>bytes_per_second;     // Nombre d'octets par seconde de musique
////  cout<<"Wave::read: Le nombre d'octets par seconde est : "<<bytes_per_second<<endl;
////  ifs>>bytes_per_sample;     // Nombre d'octets par ?chantillon
////  cout<<"Wave::read: Le nombre d'octets par ?chantillon est : "<<bytes_per_sample<<endl;
////  ifs>>depth;                // Nombre de bits par donn?e (8 ou 16)
////  cout<<"Wave::read: Le nombre de bits par donn?e est : "<<depth<<endl;
////*/
////  //chunk donn?es
////  for (i=0; i<4; i++) {
////    ifs>>data_id[i];       // Constante "data"
////    cout<<data_id[i];
////  }
////  if (!((data_id[0] == 'd') && (data_id[1] == 'a') && (data_id[2] == 't') && (data_id[3] == 'a'))) {
////    cout<<"Wave::read: Erreur, le fichier "<<fileName<<" est corrompu (pas de chaine 'data' dans le header)\n";
////    exit(-1);    
////  }
////  ifs>>data_size;       // nombre d'octet restant (i.e = taille du fichier moins 44 octets)
////
////  //LECTURE DES DONNEES
////  switch (depth) {
////    case 8:
////      samples_nb = data_size;
////      data       = new short[samples_nb];
////      char c;
////      for (i=0; i<samples_nb; i++) {
////        ifs>>c; //Lecture de 8 bits
////        data[i] = (short)c;
////      }  
////      break;
////    case 16:
////      samples_nb = data_size/2;
////      data       = new short[samples_nb];
////      for (i=0; i<samples_nb; i++) {
////        ifs>>data[i]; //Lecture de 16 bits
////      }  
////      break;
////    default:
////      cout<<"Wave::read: Erreur, la profondeur (depth = "<<depth<<") est inconnue\n";
////      exit(-1);    
////  }      
////
////  ifs.close();
//  fclose(fd);
//}    

void Wave::checkTypesSize() {
  if (sizeof(char) != 1) {
    cout<<"Wave::checkTypesSize(): Erreur : le type char n'est pas de taille 8 bits\n";
    exit(-1);
  }
  if (sizeof(short) != 2) {
    cout<<"Wave::checkTypesSize(): Erreur : le type short n'est pas de taille 16 bits\n";
    cout<<"Wave::checkTypesSize(): Erreur : le type short est pas de taille "<< 8*sizeof(int)<<" bits\n";
    exit(-1);
  }
  //cout<<"Le type int est de taille "<<sizeof(int)<<endl;
  if (sizeof(int) != 4) {
    cout<<"Wave::checkTypesSize(): Erreur : le type int n'est pas de taille 32 bits\n";
    cout<<"Wave::checkTypesSize(): Erreur : le type int est pas de taille "<< 8*sizeof(long int)<<" bits\n";
    exit(-1);
  }
}

int Wave::getFreqE() {
    return sampling_freq;
}

void Wave::filtre_passe_haut( double frequenceC ){
    int freqE = sampling_freq;
    std::cout << freqE << std::endl;
    int size = data_size;
    std::cout << size << std::endl;
    int frequenceCN = (frequenceC/(double)freqE)*size;
    std::cout << frequenceCN << std::endl;
    auto * data = new double [size];
    for (int i = 0; i < size; ++i) {
        data[i] = Char2Double(data8[i]);
    }
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(data,reel,img,size);
    writeDataFile("data_no_fph",reel,img,size);
    delete[] data;
    for (int i = 0; i < frequenceCN; ++i) {
        reel[i] = 0;
        img[i] = 0;
    }
    for (int i = size-frequenceCN; i < size; ++i) {
        reel[i] = 0;
        img[i] = 0;
    }
    writeDataFile("data_fph",reel,img,size);
    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reel,img,size);
    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {
        data_charIDFT[i] = Double2Char(dataIDFT[i]);
    }
    modifData8(data_charIDFT);
}

void Wave::filtre_passe_bas( double frequenceC ){
    int freqE = sampling_freq;
    std::cout << freqE << std::endl;
    int size = data_size;
    std::cout << size << std::endl;
    int frequenceCN = (frequenceC/(double)freqE)*size;
    std::cout << frequenceCN << std::endl;
    auto * data = new double [size];
    for (int i = 0; i < size; ++i) {
        data[i] = Char2Double(data8[i]);
    }
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(data,reel,img,size);
    writeDataFile("data_no_fph",reel,img,size);
    delete[] data;
    for (int i = frequenceCN; i < size-frequenceCN; ++i) {
        reel[i] = 0;
        img[i] = 0;
    }
    writeDataFile("data_fph",reel,img,size);
    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reel,img,size);
    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {
        data_charIDFT[i] = Double2Char(dataIDFT[i]);
    }
    modifData8(data_charIDFT);
}

int main() {
    int DO1 = 262;
    int RE = 294;
    int MI = 330;
    int FA = 349;
    int SOL = 392;
    int LA = 440;
    int SI = 494;
    int DO2 = 523;

    /*
    int nbNotes = 7;
    int secNote = 1;
    double freqE = 44100;
    int nbEchantillon = (int)secNote * (int)freqE;
    int tailleTotal = nbEchantillon*nbNotes;
    std::vector<double> notes;
    notes.emplace_back(LA);
    notes.emplace_back(SI);
    notes.emplace_back(DO2);
    notes.emplace_back(RE);
    notes.emplace_back(MI);
    notes.emplace_back(FA);
    notes.emplace_back(SOL);
    unsigned char * data_char_full = gammeChromatique(secNote,freqE,notes);
    Wave a = Wave(data_char_full,tailleTotal,1,44100);
    a.write((char*)"GammeChromatiqueLA.wav");
     */

    // LA 3 secondes
    /*
    int nbNotes = 1;
    int secNote = 6;
    double freqE = 44100;
    int nbEchantillon = (int)secNote * (int)freqE;
    int tailleTotal = nbEchantillon*nbNotes;
    std::vector<double> notes;
    notes.emplace_back(LA);
    unsigned char * data_char_full = gammeChromatique(secNote,freqE,notes);
    Wave a = Wave(data_char_full,tailleTotal,1,44100);
    a.write((char*)"La.wav");
    auto * data = new double [tailleTotal];
    for (int i = 0; i < tailleTotal; ++i) {
        data[i] = Char2Double(data_char_full[i]);
    }
    int size = 44100;
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(data,reel,img,size);
    delete[] data;
    writeDataFile("LA_data",reel,img,size);
    delete[] reel;
    delete[] img;
     */

    /*
    Wave a;
    int size = 0;
    a.read((char *)"GammePiano.wav");
    unsigned char * data_char_full;
    a.getData8(&data_char_full,&size);
    int freqE = a.getFreqE();
    auto * data = new double [size];
    for (int i = 0; i < size; ++i) {
        data[i] = Char2Double(data_char_full[i]);
    }
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(data,reel,img,size);
    delete[] data;
    writeDataFile("PIANO_data",reel,img,size);
    //reel = new double [size];
    //img = new double [size];
    //readDataFile("PIANO_data",reel,img);
    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reel,img,size);
    delete[] reel;
    delete[] img;
    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {
        data_charIDFT[i] = Double2Char(dataIDFT[i]);
    }
    Wave b = Wave(data_charIDFT,size,1,freqE);
    b.write((char*)"GaPiano.wav");
    delete[] dataIDFT;
    delete[] data_charIDFT;
    */

    /*
    Wave a;
    a.read((char *)"BrokenGlass.wav");
    unsigned char * data;
    int size;
    a.getData8(&data,&size);
    auto * signal = new double [size];
    for (int i = 0; i < size; ++i) {signal[i] = Char2Double(data[i]);}
    auto * signalFilter = new double [size];
    filtreButterwoth(signal,signalFilter,size,a.getFreqE(),4000);
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(signalFilter,reel,img,size);
    writeDataFile("ButterPop",reel,img,size);
    for (int i = 0; i < size; ++i) {data[i] = Double2Char(signalFilter[i]);}
    a.modifData8(data);
    a.write((char *)"BrokenGlass_ButterWorth.wav");
    */


    Wave a;
    a.read((char *)"Full_Note.wav");
    unsigned char * data;
    int size;
    a.getData8(&data,&size);
    auto * signal = new double [size];
    for (int i = 0; i < size; ++i) {signal[i] = Char2Double(data[i]);}
    auto * reel = new double [size];
    auto * img = new double [size];
    DFT(signal,reel,img,size);
    //writeDataFile("Full_Note_DATA",reel,img,size);
    auto * signalFilter = new double [size];
    filtreButterwoth(signal,signalFilter,size,a.getFreqE(),DO1+1);
    //filtreButterwoth(signal,signalFilter,size,2);
    //filtreBessel(signal,signalFilter,size);
    reel = new double [size];
    img = new double [size];
    DFT(signalFilter,reel,img,size);
    writeDataFile("Full_Note_DATA_Filtre_Butterwoth",reel,img,size);
    for (int i = 0; i < size; ++i) {data[i] = Double2Char(signalFilter[i]);}
    a.modifData8(data);
    //a.write((char *)"LADOMI_ButterWorth.wav");
    a.write((char *)"Full_Note_Butterwoth.wav");


    //a.filtre_passe_haut(5000);
    //a.write((char *)"BrokenGlass_FPH.wav");
    //a.filtre_passe_bas(6000);
    //a.write((char *)"BrokenGlass_FPB.wav");


    /*
    int secNote = 5;
    double freqE = 8000;
    int size = secNote * (int)freqE;
    auto * reelDOMISOL = new double [size];
    auto * imgDOMISOL = new double [size];

    std::vector<double> notes;
    notes.emplace_back(DO1);
    notes.emplace_back(RE);
    notes.emplace_back(MI);
    notes.emplace_back(SOL);
    notes.emplace_back(FA);
    notes.emplace_back(LA);
    notes.emplace_back(SI);
    notes.emplace_back(DO2);
    accord(notes,secNote,reelDOMISOL,imgDOMISOL,freqE);

    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reelDOMISOL,imgDOMISOL,size);
    delete[] reelDOMISOL;
    delete[] imgDOMISOL;

    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {data_charIDFT[i] = Double2Char(dataIDFT[i]);}
    delete[] dataIDFT;
    Wave b = Wave(data_charIDFT,size,1,(int)freqE);
    delete[] data_charIDFT;
    b.write((char*)"Full_Note.wav");
    */



    /*
    int secNote = 5;
    double freqE = 8000;
    int size = secNote * (int)freqE;
    auto * reelDOMISOL = new double [size];
    auto * imgDOMISOL = new double [size];

    std::vector<double> notes;
    notes.emplace_back(LA);
    notes.emplace_back(DO1);
    notes.emplace_back(MI);
    accord(notes,secNote,reelDOMISOL,imgDOMISOL,freqE);

    writeDataFile("accord_data",reelDOMISOL,imgDOMISOL,size);

    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reelDOMISOL,imgDOMISOL,size);
    delete[] reelDOMISOL;
    delete[] imgDOMISOL;

    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {data_charIDFT[i] = Double2Char(dataIDFT[i]);}
    delete[] dataIDFT;
    Wave b = Wave(data_charIDFT,size,1,(int)freqE);
    //b.write((char*)"LADOMI.wav");
    Wave * output = b.filtre_passe_haut(450);
    output->write((char*)"LADOMI_PH.wav");
    */

    /*
    auto * reel = new double [size];
    auto * img = new double [size];
    readDataFile("accord_data",reel,img);
    auto * dataIDFTTrans = new double[size];
    IDFT(dataIDFTTrans,reel,img,size);
    auto * data_charIDFTTrans = new unsigned char [size];
    for (int i = 0; i < size; ++i) {data_charIDFTTrans[i] = Double2Char(dataIDFTTrans[i]);}
    delete[] dataIDFTTrans;
    Wave c = Wave(data_charIDFTTrans,size,1,(int)freqE);
    delete[] data_charIDFTTrans;
    //c.write((char*)"LADOMI_Trans.wav");
     */


    /*
    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reel,img,size);
    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {
        data_charIDFT[i] = Double2Char(dataIDFT[i]);
    }
    Wave b = Wave(data_charIDFT,tailleTotal,1,44100);
    b.write((char*)"NoteIDFT.wav");
     */

    /*
    int nbNotes = 1;
    double secNote = 5;
    double freqE = 44100;
    int nbEchantillon = (int)secNote * (int)freqE;
    int tailleTotal = nbEchantillon*nbNotes;
    int size = 44100;
    auto * reel = new double [size];
    auto * img = new double [size];
    readDataFile("DO_MI_SOL",reel,img,size);
    auto * dataIDFT = new double[size];
    IDFT(dataIDFT,reel,img,size);
    delete[] reel;
    delete[] img;
    auto * data_charIDFT = new unsigned char [size];
    for (int i = 0; i < size; ++i) {
        data_charIDFT[i] = Double2Char(dataIDFT[i]);
    }
    delete[] dataIDFT;
    Wave b = Wave(data_charIDFT,tailleTotal,1,size);
    delete[] data_charIDFT;
    b.write((char*)"NoteIDFT.wav");
     */
    return 0;
}