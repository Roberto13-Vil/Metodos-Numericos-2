//newton
#include <stdio.h>
#include <math.h>

double m_cofactor[10][10];
double m_transpuesta[10][10];
double jacobinainversa1[10][10];
double jacobina[10][10];
double evalFunc[10][10];
double multJPorF[10][10];
double deter=0.0;

double absolute(double x);
void multiplicacion(double x[10][10],int f1, int c1, double y[10][10],int f2, int c2);
double determinante(double matriz[10][10], int orden);
double cofactor(double matriz[10][10], int orden,int fila, int columna);
void MatrizCofactor(double matriz[10][10], int orden);
void MatrizInversa(double matriz[10][10], int orden);

// Sistemas
void primerSistema();
void segundoSistema();
void tercerSistema();
void cuartoSistema();

int main(){
    int opcion=100;
    do{
        printf("Metodo de Newton Rapson\n\n");
        printf("1. Sistema de ecuaciones:\n");
        printf("f1(x, y) = x^2+xy-10=0\n");
        printf("f2(x, y) = y+3xy^2-50=0\n\n");

        printf("2. Sistema de ecuaciones:\n");
        printf("f1(x, y) = x^2+y^2-9=0\n");
        printf("f2(x, y) = -e^x-2y-3=0\n\n");

        printf("3. Sistema de ecuaciones:\n");
        printf("f1(x, y,z) = 2x^2-4x+y^2+3z^2+6z+2=0\n");
        printf("f2(x, y,z) = x^2+y^2-2y+2z^2-5=0\n");
        printf("f2(x, y,z) = 3x^2-12x+y^2-3z^2+8=0\n\n");

        printf("4. Sistema de ecuaciones:\n");
        printf("f1(x, y,z) = x^2-4x+y^2=0\n");
        printf("f2(x, y,z) =x^2-x-12y+1=0 \n");
        printf("f2(x, y,z) =3x^2-12x+y^2-3z^2+8=0\n\n");

        printf("5. Salir\n\n");

        printf("Selecciona una opcion: ");
        scanf("%d", &opcion);

        system("cls");
        switch(opcion){
        case 1:
            primerSistema();
            break;
        case 2:
            segundoSistema();
            break;
        case 3:
            tercerSistema();
            break;
        case 4:
            cuartoSistema();
            break;
        case 5:
            opcion=5;
            break;
        default:
            printf("Esa opción no es valida intente de nuevo.");
            system("pause");
        }
        system("cls");
    }while(opcion!=5);
    return 0;
}

void primerSistema(){
    char mas='s';
    char aux='s';
    do{
    printf("El sistema de ecuaciones es:\n");
    printf("f1(x, y) = x^2+xy-10=0\n");
    printf("f2(x, y) = y+3xy^2-50=0\n\n");

    // derivada parcial de la primera funcion
    printf("f1'x=2x+y\t");
    printf("f1'y=x\t\n");

    // derivada parcial de la segunda funcion
    printf("f2'x=3y^2\t");
    printf("f2'y=6yx+1\t\n");

    double puntos[2];
    double puntos2[2];
    double errorauxpuntos[2];
    double maxaux1;
    double maxaux2;
    double error;
    double tol;
    int it;

    printf("Ingresa el primer punto: ");
    scanf("%lf", &puntos[0]);
    printf("\nIngresa el segundo punto: ");
    scanf("%lf", &puntos[1]);
    printf("\nIngresa las iteraciones: ");
    scanf("%d", &it);
    printf("Ingresa la tolerancia: ");
    scanf("%lf", &tol);

    for(int k=0; k<it; k++){
    jacobina[0][0] = 2*puntos[0]+puntos[1];
    jacobina[0][1] = puntos[0];
    jacobina[1][0] = 3*(puntos[1]*puntos[1]);
    jacobina[1][1] = 6*puntos[1]*puntos[0]+1;

    evalFunc[0][0] = (puntos[0]*puntos[0])+(puntos[0]*puntos[1])-10;
    evalFunc[1][0] = puntos[1]+(3*puntos[0]*(puntos[1]*puntos[1]))-50;

    //inversa
    deter = determinante(jacobina, 2);
    MatrizCofactor(jacobina, 2);
    //printf("Matriz inversa:\n");
    MatrizInversa(m_cofactor, 2);

    //printf("Multiplicacion: \n");
    multiplicacion(jacobinainversa1, 2,2,evalFunc,2,1);

    //printf("%f ", multJPorF[0][0]);
    //printf("%f", multJPorF[1][0]);

    puntos2[0] = puntos[0] - multJPorF[0][0];
    puntos2[1] = puntos[1] - multJPorF[1][0];

    //printf("\npuntos nuevos: \n");
    for(int i=0; i<2; i++){
        errorauxpuntos[i] = absolute(puntos2[i]-puntos[i]);
    }
    //printf("\nError: ");

    if(errorauxpuntos[0]>errorauxpuntos[1]){
        maxaux1 = errorauxpuntos[0];
    }
    else{
        maxaux1 = errorauxpuntos[1];
    }
    if(absolute(puntos2[0])>absolute(puntos2[1])){
        maxaux2 = puntos2[0];
    }
    else{
        maxaux2 = puntos2[1];
    }
    // .101369
    error = maxaux1/absolute(maxaux2);

    puntos[0] = puntos2[0];
    puntos[1] = puntos2[1];
    if(error<tol){
        printf("Fin por toleracia: \n");
        printf("Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f", puntos2[0], puntos2[1]);
        break;
    }
    else if(it==k+1){
        printf("Fin por Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f", puntos2[0], puntos2[1]);
        break;
    }
    printf("\n\n");
    printf("ITERACION: %d\n", k+1);
    printf("ERROR: %.8f\n", error);
    printf("PUNTOS: x: %f y: %f\n", puntos2[0], puntos2[1]);
    }
    printf("\n¿Quieres ingresar nuevos datos? s o n\n");
    fflush(stdin);
    scanf(" %c", &aux);

    //printf("%c", aux);
    if(aux=='n'){
        mas='n';
    }
    system("cls");
    }while(mas=='s');
}

void segundoSistema(){
    char mas='s';
    char aux='s';
    do{
    printf("El sistema de ecuaciones es:\n");
    printf("f1(x, y) = x^2+y^2-9=0\n");
    printf("f2(x, y) = -e^x-2y-3=0\n\n");

    // derivada parcial de la primera funcion
    printf("f1'x=2x \t");
    printf("f1'y=2y\t\n");

    // derivada parcial de la segunda funcion
    printf("f2'x=-e^x\t");
    printf("f2'y=-2\t\n");

    double puntos[2];
    double puntos2[2];
    double errorauxpuntos[2];
    double maxaux1;
    double maxaux2;
    double error;
    double tol;
    int it;

    printf("Ingresa el primer punto: ");
    scanf("%lf", &puntos[0]);
    printf("\nIngresa el segundo punto: ");
    scanf("%lf", &puntos[1]);
    printf("\nIngresa las iteraciones: ");
    scanf("%d", &it);
    printf("Ingresa la tolerancia: ");
    scanf("%lf", &tol);

    for(int k=0; k<it; k++){
    jacobina[0][0] = 2*puntos[0];
    jacobina[0][1] = 2*puntos[1];
    jacobina[1][0] = -exp(puntos[0]);
    jacobina[1][1] = -2;

    evalFunc[0][0] = (puntos[0]*puntos[0])+(puntos[1]*puntos[1])-9;
    evalFunc[1][0] = -exp(puntos[0])-(2*puntos[1])-3;

    //inversa
    deter = determinante(jacobina, 2);
    MatrizCofactor(jacobina, 2);
    //printf("Matriz inversa:\n");
    MatrizInversa(m_cofactor, 2);

    //printf("Multiplicacion: \n");
    multiplicacion(jacobinainversa1, 2,2,evalFunc,2,1);

    //printf("%f ", multJPorF[0][0]);
    //printf("%f", multJPorF[1][0]);

    puntos2[0] = puntos[0] - multJPorF[0][0];
    puntos2[1] = puntos[1] - multJPorF[1][0];

    //printf("\npuntos nuevos: \n");
    for(int i=0; i<2; i++){
        errorauxpuntos[i] = absolute(puntos2[i]-puntos[i]);
    }
    //printf("\nError: ");

    if(errorauxpuntos[0]>errorauxpuntos[1]){
        maxaux1 = errorauxpuntos[0];
    }
    else{
        maxaux1 = errorauxpuntos[1];
    }
    if(absolute(puntos2[0])>absolute(puntos2[1])){
        maxaux2 = puntos2[0];
    }
    else{
        maxaux2 = puntos2[1];
    }
    // .101369
    error = maxaux1/absolute(maxaux2);

    puntos[0] = puntos2[0];
    puntos[1] = puntos2[1];
    if(error<tol){
        printf("Fin por toleracia: \n");
        printf("Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f", puntos2[0], puntos2[1]);
        break;
    }
    else if(it==k+1){
        printf("Fin por Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f", puntos2[0], puntos2[1]);
        break;
    }
    printf("\n\n");
    printf("ITERACION: %d\n", k+1);
    printf("ERROR: %.8f\n", error);
    printf("PUNTOS: x: %f y: %f\n", puntos2[0], puntos2[1]);
    }
    printf("\n¿Quieres ingresar nuevos datos? s o n\n");
    fflush(stdin);
    scanf(" %c", &aux);

    //printf("%c", aux);
    if(aux=='n'){
        mas='n';
    }
    system("cls");
    }while(mas=='s');
}

void tercerSistema(){
    char mas='s';
    char aux='s';
    do{
    printf("El sistema de ecuaciones es:\n");
    printf("f1(x, y,z) = 2x^2-4x+y^2+3z^2+6z+2=0\n");
    printf("f2(x, y,z) = x^2+y^2-2y+2z^2-5=0\n");
    printf("f2(x, y,z) = 3x^2-12x+y^2-3z^2+8=0\n\n");

    // derivada parcial de la primera funcion
    printf("f1'x=4x-4\t");
    printf("f1'y=2y \t");
    printf("f1'z=6z+6\t\n");

    // derivada parcial de la segunda funcion
    printf("f2'x=2x \t");
    printf("f2'y=2y-2\t");
    printf("f2'z=4z\t\n");

    // derivada parcial de la tercera funcion
    printf("f3'x=6x-12\t");
    printf("f3'y=2y \t");
    printf("f3'z=-6z\t\n");

    double puntos[3];
    double puntos2[3];
    double errorauxpuntos[3];
    double maxaux1;
    double maxaux2;
    double error;
    double tol;
    int it;

    printf("Ingresa el primer punto: ");
    scanf("%lf", &puntos[0]);
    printf("\nIngresa el segundo punto: ");
    scanf("%lf", &puntos[1]);
    printf("\nIngresa el tercer punto: ");
    scanf("%lf", &puntos[2]);
    printf("\nIngresa las iteraciones: ");
    scanf("%d", &it);
    printf("Ingresa la tolerancia: ");
    scanf("%lf", &tol);

    for(int k=0; k<it; k++){
    jacobina[0][0] = 4*puntos[0]-4;
    jacobina[0][1] = 2*puntos[1];
    jacobina[0][2] = 6*puntos[2]+6;
    jacobina[1][0] = 2*puntos[0];
    jacobina[1][1] = 2*puntos[1]-2;
    jacobina[1][2] = 4*puntos[2];
    jacobina[2][0] = 6*puntos[0]-12;
    jacobina[2][1] = 2*puntos[1];
    jacobina[2][2] = -6*puntos[2];

    evalFunc[0][0] =2*pow(puntos[0],2)-4*puntos[0]+pow(puntos[1],2)+3*pow(puntos[2],2)+6*puntos[2]+2;
    evalFunc[1][0] =pow(puntos[0],2)+pow(puntos[1],2)-2*puntos[1]+2*pow(puntos[2],2)-5;
    evalFunc[2][0] =3*pow(puntos[0],2)-12*puntos[0]+pow(puntos[1],2)-3*pow(puntos[2],2)+8;

    //inversa
    deter = determinante(jacobina, 3);
    MatrizCofactor(jacobina, 3);
    //printf("Matriz inversa:\n");
    MatrizInversa(m_cofactor, 3);

    //printf("Multiplicacion: \n");
    multiplicacion(jacobinainversa1, 3,3,evalFunc,3,1);

    //printf("%f ", multJPorF[0][0]);
    //printf("%f", multJPorF[1][0]);

    puntos2[0] = puntos[0] - multJPorF[0][0];
    puntos2[1] = puntos[1] - multJPorF[1][0];
    puntos2[2] = puntos[2] - multJPorF[2][0];

    //printf("\npuntos nuevos: \n");
    for(int i=0; i<3; i++){
        errorauxpuntos[i] = absolute(puntos2[i]-puntos[i]);
    }
    //printf("\nError: ");

    if(errorauxpuntos[0]>errorauxpuntos[1] && errorauxpuntos[0]>errorauxpuntos[2]){
        maxaux1 = errorauxpuntos[0];
    }
    else if(errorauxpuntos[1]>errorauxpuntos[0] && errorauxpuntos[1]>errorauxpuntos[2]){
        maxaux1 = errorauxpuntos[1];
    }
    else{
        maxaux1 = errorauxpuntos[2];
    }

    if(absolute(puntos2[0])>absolute(puntos2[1])&&absolute(puntos2[0])>absolute(puntos2[2])){
        maxaux2 = puntos2[0];
    }
    else if(absolute(puntos2[1])>absolute(puntos2[0])&&absolute(puntos2[1])>absolute(puntos2[2])){
        maxaux2 = puntos2[1];
    }
    else{
        maxaux2 = puntos2[2];
    }
    // .101369
    error = maxaux1/absolute(maxaux2);

    puntos[0] = puntos2[0];
    puntos[1] = puntos2[1];
    puntos[2] = puntos2[2];

    if(error<tol){
        printf("Fin por toleracia: \n");
        printf("Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f z: %f", puntos2[0], puntos2[1], puntos2[2]);
        break;
    }
    else if(it==k+1){
        printf("Fin por Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f z: %f", puntos2[0], puntos2[1], puntos2[2]);
        break;
    }
    printf("\n\n");
    printf("ITERACION: %d\n", k+1);
    printf("ERROR: %.8f\n", error);
    printf("PUNTOS: x: %f y: %f z: %f\n", puntos2[0], puntos2[1], puntos2[2]);
    }
    printf("\n¿Quieres ingresar nuevos datos? s o n\n");
    fflush(stdin);
    scanf(" %c", &aux);

    //printf("%c", aux);
    if(aux=='n'){
        mas='n';
    }
    system("cls");
    }while(mas=='s');
}

void cuartoSistema(){
    char mas='s';
    char aux='s';
    do{
    printf("El sistema de ecuaciones es:\n");
    printf("f1(x, y,z) = x^2-4x+y^2=0\n");
    printf("f2(x, y,z) =x^2-x-12y+1=0 \n");
    printf("f2(x, y,z) =3x^2-12x+y^2-3z^2+8=0\n");

    // derivada parcial de la primera funcion
    printf("f1'x=2x-4\t");
    printf("f1'y=2y\t");
    printf("f1'z=0\t\n");

    // derivada parcial de la segunda funcion
    printf("f2'x=2x-1 \t");
    printf("f2'y=-12\t");
    printf("f2'z=0\t\n");

    // derivada parcial de la tercera funcion
    printf("f3'x=6x-12\t");
    printf("f3'y=2y \t");
    printf("f3'z=-6z\t\n");

    double puntos[3];
    double puntos2[3];
    double errorauxpuntos[3];
    double maxaux1;
    double maxaux2;
    double error;
    double tol;
    int it;

    printf("Ingresa el primer punto: ");
    scanf("%lf", &puntos[0]);
    printf("\nIngresa el segundo punto: ");
    scanf("%lf", &puntos[1]);
    printf("\nIngresa el tercer punto: ");
    scanf("%lf", &puntos[2]);
    printf("\nIngresa las iteraciones: ");
    scanf("%d", &it);
    printf("Ingresa la tolerancia: ");
    scanf("%lf", &tol);

    for(int k=0; k<it; k++){
    jacobina[0][0] = 2*puntos[0]-4;
    jacobina[0][1] = 2*puntos[1];
    jacobina[0][2] = 0;
    jacobina[1][0] = 2*puntos[0]-1;
    jacobina[1][1] = -12;
    jacobina[1][2] = 0;
    jacobina[2][0] = 6*puntos[0]-12;
    jacobina[2][1] = 2*puntos[1];
    jacobina[2][2] = -6*puntos[2];

    evalFunc[0][0] = pow(puntos[0],2)-4*puntos[0]+pow(puntos[1],2);
    evalFunc[1][0] = pow(puntos[0],2)-puntos[0]-12*puntos[1]+1;
    evalFunc[2][0] = 3*pow(puntos[0],2)-12*puntos[0]+pow(puntos[1],2)-3*pow(puntos[2],2)+8;

    //inversa
    deter = determinante(jacobina, 3);
    MatrizCofactor(jacobina, 3);
    //printf("Matriz inversa:\n");
    MatrizInversa(m_cofactor, 3);

    //printf("Multiplicacion: \n");
    multiplicacion(jacobinainversa1, 3,3,evalFunc,3,1);

    //printf("%f ", multJPorF[0][0]);
    //printf("%f", multJPorF[1][0]);

    puntos2[0] = puntos[0] - multJPorF[0][0];
    puntos2[1] = puntos[1] - multJPorF[1][0];
    puntos2[2] = puntos[2] - multJPorF[2][0];

    //printf("\npuntos nuevos: \n");
    for(int i=0; i<3; i++){
        errorauxpuntos[i] = absolute(puntos2[i]-puntos[i]);
    }
    //printf("\nError: ");

    if(errorauxpuntos[0]>errorauxpuntos[1] && errorauxpuntos[0]>errorauxpuntos[2]){
        maxaux1 = errorauxpuntos[0];
    }
    else if(errorauxpuntos[1]>errorauxpuntos[0] && errorauxpuntos[1]>errorauxpuntos[2]){
        maxaux1 = errorauxpuntos[1];
    }
    else{
        maxaux1 = errorauxpuntos[2];
    }

    if(absolute(puntos2[0])>absolute(puntos2[1])&&absolute(puntos2[0])>absolute(puntos2[2])){
        maxaux2 = puntos2[0];
    }
    else if(absolute(puntos2[1])>absolute(puntos2[0])&&absolute(puntos2[1])>absolute(puntos2[2])){
        maxaux2 = puntos2[1];
    }
    else{
        maxaux2 = puntos2[2];
    }
    // .101369
    error = maxaux1/absolute(maxaux2);

    puntos[0] = puntos2[0];
    puntos[1] = puntos2[1];
    puntos[2] = puntos2[2];

    if(error<tol){
        printf("Fin por toleracia: \n");
        printf("Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f z: %f", puntos2[0], puntos2[1], puntos2[2]);
        break;
    }
    else if(it==k+1){
        printf("Fin por Iteraciones: %d\n", k+1);
        printf("error: %.8f\n", error);
        printf("puntos: x: %f y: %f z: %f", puntos2[0], puntos2[1], puntos2[2]);
        break;
    }
    printf("\n\n");
    printf("ITERACION: %d\n", k+1);
    printf("ERROR: %.8f\n", error);
    printf("PUNTOS: x: %f y: %f z: %f\n", puntos2[0], puntos2[1], puntos2[2]);
    }
    printf("\n¿Quieres ingresar nuevos datos? s o n\n");
    fflush(stdin);
    scanf(" %c", &aux);

    //printf("%c", aux);
    if(aux=='n'){
        mas='n';
    }
    system("cls");
    }while(mas=='s');
}

double absolute(double x){
    if(x<0){
        return x*-1;
    }
    else
        return x;
}
double determinante(double matriz[10][10], int orden){
   double det = 0.0;

   if (orden == 1) {
      det = matriz[0][0];
   } else {
      for (int j = 0; j < orden; j++) {
         det = det + matriz[0][j] * cofactor(matriz, orden, 0, j);
      }
   }

   return det;
}

double cofactor(double matriz[10][10], int orden,int fila, int columna){
   double submatriz[10][10];
   int n = orden - 1;
   int i, j;

   int x = 0;
   int y = 0;
   for (i = 0; i < orden; i++) {
      for (j = 0; j < orden; j++) {
         if (i != fila && j != columna) {
            submatriz[x][y] = matriz[i][j];
            y++;
            if (y >= n) {
               x++;
               y = 0;
            }
         }
      }
   }
   return pow(-1.0, fila + columna) * determinante(submatriz, n);
}

void multiplicacion(double x[10][10],int f1, int c1, double y[10][10],int f2, int c2){
    for(int i=0; i<c2; i++){
        for(int j=0; j<f1; j++){
            double suma=0;
            for(int z=0; z<c1; z++){
                suma+= x[j][z]*y[z][i];
            }
        multJPorF[j][i] = suma;
        }
    }
}

void MatrizCofactor(double matriz[10][10], int orden){
   for(int i=0; i<orden; i++) {
       for(int j=0; j<orden; j++){
           m_cofactor[i][j]= cofactor(matriz, orden, i, j);
       }
    }
}

void MatrizInversa(double matriz[10][10], int orden){
    for(int i=0; i<orden; i++) {
       for(int j=0; j<orden; j++){
           m_transpuesta[j][i] = m_cofactor[i][j];
       }
    }
    for(int i=0; i<orden; i++) {
       for(int j=0; j<orden; j++){
           jacobinainversa1[i][j] = m_transpuesta[i][j]/deter;
       }
    }
}

