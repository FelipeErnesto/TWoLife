#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <R.h>
#include <Rmath.h>
#include "paisagem.cpp" 
#include "individuo.cpp" 
#include <sstream>
#include <string>

using namespace std;

/** \file TWoLife.cpp \brief Arquivo usado na integração R/C
 * 
 * Este arquivo define a função TWoLife em C, que será chamada a partir do R */


/** Funcao que recebe os parâmetros da simulação e a executa,
 // baseada no texto em http://users.stat.umn.edu/~geyer/rc/
 //
 // NOTE que funções em C++ precisam ser indicadas como extern "C"
 // para poderem ser acessadas facilmente pela interface do R!
 //
 // Veja a descrição da classe \ref paisagem para o significado dos parâmetros
 // */
extern "C" void TWoLife (double * raio, int * N, double * angulo_visada, double * passo, double * move,
						 double * taxa_basal, double * taxa_morte, double * incl_b, double * incl_d,
						 int * numb_cells, double * cell_size, int * land_shape, int * density_type, 
						 double * death_mat, int * inipos, int * bound_condition, int * scape, double * tempo, int * nPop, double * x, double * y, int * outCode)
{
	// This sequence creates an attribute containing the output file name. The template is output-00000.txt.
	string fileNAME = "output-00000.txt";
	stringstream tmps;
	tmps<<outCode[0];
	string addToName = tmps.str();
	int fnSize = fileNAME.size();
    int tmpsSize = addToName.size();
	fileNAME.erase(fileNAME.begin()+fnSize-4-tmpsSize,fileNAME.begin()+fnSize-4);
	fileNAME.insert(fnSize-4-tmpsSize,addToName);
	
	GetRNGstate(); /* (mudar para doxygen):  este comando chama o engine de numeros aleatorios do R
					Por causa dela nossa biblioteca nao eh standalone */
	
	paisagem* floresta = new paisagem(raio[0], N[0], angulo_visada[0], passo[0], 
									  move[0], taxa_basal[0], taxa_morte[0], incl_b[0], 
									  incl_d[0], numb_cells[0], cell_size[0], land_shape[0],
									  density_type[0], death_mat[0], inipos[0], bound_condition[0], 
									  scape);
	
	ofstream outputSIM; // ofstream for the output file
	outputSIM.open(fileNAME.c_str());
	for(unsigned int i=0; i<floresta->conta_individuos();i++)
	{
		outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
	}
	
	while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
	{
		int ind_neo = floresta->sorteia_individuo();
		int acao = floresta->sorteia_acao(ind_neo);
		individuo* atestado = new individuo(*floresta->get_individuos(ind_neo));
		individuo::subtract_id();
		floresta->realiza_acao(acao, ind_neo);
		double t_ant = floresta->tempo_do_mundo;
		floresta->update(acao, ind_neo, atestado);
		if(t_ant < (int)floresta->tempo_do_mundo)
		{
			for(unsigned int i=0; i<floresta->conta_individuos();i++)
			{
				outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
			}
		}
	}
	if(floresta->conta_individuos()==0){outputSIM << floresta->tempo_do_mundo << " " << "NA" << " " << "NA" << " " << "NA" << endl;}
	outputSIM.close(); //end of output file
	
	ofstream arqm;
	ofstream arqe;
	arqm.open ("migracao.txt");
	arqe.open ("extincao.txt");

	int numb_patches = floresta->get_numb_patches();
	arqm << numb_patches<<endl;
	arqe << numb_patches<<endl;
	for (int j = 1; j<numb_patches+1; j++)
	{
		arqm << floresta->get_migracao(j)<<",";
		arqe << floresta->get_extincao(j)<<",";
	}
	arqm.close();
	arqe.close();
	
	*nPop = floresta->conta_individuos();
	for (int i =0; i < *nPop; i ++) {
		x[i] = floresta->get_individuos(i)->get_x();
		y[i] = floresta->get_individuos(i)->get_y();
	} //DUVIDA: porque x[i] e y[i] nao tem asterisco antes?
	delete floresta;
	PutRNGstate();
}


