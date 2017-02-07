#ifndef PAISAGEM_H
#define PAISAGEM_H
#include "individuo.hpp"
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<cmath>
#include <forward_list>
///** Dimensão maxima do mapa, usada para rotinas de fragmentação TBI */
#define dim 10000 /*Em algum momento pode ser alterado para um argumento do
		    construtor. No momento não é prioritário. */

using namespace std;

/** \brief A classe paisagem implementa o mundo onde os agentes estão.
 *
 * Esta classe é responsável pela criação dos indivíduos, por manter o relógio do mundo, e por estabelecer a comunicação entre os indivíduos \sa \ref paisagem::update */
class paisagem
{
private:

    //propriedades privadas
    const double tamanho; //tamanho do lado do mundo
	const unsigned long N; // numero de inds no comeco da simulacao
	vector <individuo*> popIndividuos;//vetor dos individuos da populacao
    const int numb_cells; // em quantos pedacos cada dimensao do mundo esta dividida (para fragmentacao)
    const double cell_size;// tamanho do lada de uma célula ou tamanho do pixel da imagem usada para a classificação da paisagem
	const int landscape_shape; // paramâmetro do modelo relacionado à forma da paisagem (por enquanto, 1 = "quadrada" ou  0 = "circular")
	const int boundary_condition; // tipo de condição de contorno (0 = absortiva, 1 = periódica, 2 = reflexiva)
	int landscape[dim][dim];//[linha][coluna] temporariamente substituido or valor fixo
	int patches[dim][dim]; //Matriz com a determinação do fragmento a que cada pixel pertence. 0 para matriz; i>=1 para fragmentos
	forward_list<forward_list<double> > dmatrix;//matriz de distancias
	int numb_patches; //Número de fragmentos encontrados. Desconsidera-se a matriz.
	int* migracao;
	int* patch_pop;
	int* extincao;
	int* patch_area;
	const int initialPos;
	//metodos privados
	void populating(
					const double raio, 
					/** Numero de individuos que começam a simulação espalhados no ambiente */
					const int N, 
					/** Ângulo de visada dos individuos */
					const double angulo_visada,
					/** Passo de caminhada dos individuos */
					const double passo, 
					/** Taxa de movimentação dos individuos */
					const double move, 
					/** Taxa de nascimento de um individuo no habitat favoravel e sem vizinhos */
					const double taxa_basal, 
					/** Taxa de morte dos individuos */
					const double taxa_morte,		
					/** Inclinação da denso-dependência na natalidade */
					const double incl_b,
					/** Inclinação da denso-dependência na mortalidade */
					const double incl_d,
					/** Constante de acrescimo na taxa basal de morte na matriz em relação à essa taxa basal no habitat */ 
					const double death_m,
					/** Tipo de densidade (0 = GLOBAL, 1 = LOCAL) */
					const int dens_type
					);
					
    void atualiza_vizinhos(int i) const;//contabilizador de vizinhos
    void atualiza_habitat(individuo * const ind) const;//vai informar o individuo em que tipo de habitat ele esta
    void initialize_patch(individuo * const ind) const;//Vai informar o indivíduo recém criado em que fragmento ele está
    void atualiza_patch(individuo * const ind) const;//vai informar o individuo em que fragmento ele está
    void atualiza_migracao(individuo * const ind) const;
    void atualiza_extincao(int label) const;
    void initialize_dmatrix(); //Por algum motivo no consigo colocar const aqui :/
    //int define_tempo();
	void apply_boundary(individuo * const ind); //const; // metodo para aplicação da condicao de contorno
    		

public:
	//vector <individuo*> popIndividuos;
	
	/** Contador de quanto tempo já transcorreu no mundo simulado */
    double tempo_do_mundo;
    
	//metodos publicos
	/** Construtor da classe paisagem */
    paisagem(
			/** Raio no qual os individuos percebem vizinhos. Deve ser menor que o tamanho do mundo */
			const double raio, 
			/** Numero de individuos que começam a simulação espalhados no ambiente */
			const int N, 
			/** Ângulo de visada dos individuos */
			const double angulo_visada,
            /** Passo de caminhada dos individuos */
			const double passo, 
			/** Taxa de movimentação dos individuos */
			const double move, 
			/** Taxa de nascimento de um individuo no habitat favoravel e sem vizinhos */
			const double taxa_basal, 
			/** Taxa de morte dos individuos */
			const double taxa_morte,
			/** Inclinação da denso-dependência na natalidade */
			const double incl_b,
			/** Inclinação da denso-dependência na mortalidade */
			const double incl_d,
			/** Número de células de cada lado da grade que representa a paisagem */
			const int numb_cells,
			/** Tamanho do lado da célula (ou tamanho do pixel) */
			const double cell_size,			 
			/** Forma da paisagem */
			const int land_shape,
			/** Tipo de densidade ("g" = GLOBAL, "l" = LOCAL) */
			const int density_type,
			/** Constante de acrescimo na taxa basal de morte na matriz em relação à essa taxa basal no habitat */ 
			const double death_mat,
			/** How individuals are initially set into the landscpae **/
			const int inipos,
			/** Condição de contorno */
			const int bound_condition,
			/** Vetor de cobertura de habitat na paisagem */
			int scape[]			
			); //construtor

	/** Atualiza as variáveis de todos os indivíduos (ver individuo::set_vizinhos, individuo::set_habitat e individuo::update) e escolhe uma ação para ser executada. Executa a ação e atualiza o tempo do mundo de acordo \sa \ref Introdução */
    void update(int acao, int ind_chosen, individuo* atestado);//atualizador
    int sorteia_individuo();
    int sorteia_acao(const int lower){return this->popIndividuos[lower]->sorteia_acao();}
	void realiza_acao(int acao, int lower);//vai pegar os tempos de cada individuo e informa qual foi o escolhido e manda ele fazer
	/** Retorna o número total de indivíduos na paisagem */
    const int conta_individuos() const{return popIndividuos.size();}
	/** Retorna um vetor contendo todos os indivíduos na paisagem */
    individuo* get_individuos(int i) const {return popIndividuos[i];}
	/** Retorna o número de espécies existentes na paisagem \ref TBI */
    const int conta_especies() const;
	/** Retorna o tamanho da paisagem (definido no construtor) */
    const double get_tamanho() const {return this->tamanho;}
	
	double calcDist(const individuo* a1, const individuo* a2) const;
	
	double calcDensity(const individuo* ind1) const;

    /** Retorna false se o indivíduo estava no ambiente no passod de tempo 0, e true se ele nasceu durante a simulação.
	 * Usado para pintar os indivíduos nascidos de um cor diferente dos individuos originais */
    const bool nascido(individuo * const ind) const {return ind->get_id() > this->N;}
    
    void find_patches(int x, int y, int current_label);
	
    int get_migracao(int i){return migracao[i];}
    int get_extincao(int i){return extincao[i];}
    int get_numb_patches(){return numb_patches;}
    double get_dist(int ind1, int ind2);
};

#endif // PAISAGEM_H
