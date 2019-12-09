#include <iostream>
#include <math.h>
#include "Eigen/Eigen/Eigen"
#include "freeglut/include/GL/glut.h"

using namespace Eigen;

//parametri prozora
static int window_width = 900;
static int window_height = 600;
static int window_x = 250;
static int window_y = 100;

//Ojlerovi uglovi za prvo stanje koordiatnog sistema i drugo
static double alfa_1 = M_PI/2, beta_1 = 3*M_PI / 4, gamma_1 = M_PI / 2;
static double alfa_2 = -M_PI/4, beta_2 = 2 * M_PI / 3, gamma_2 = -4*M_PI/3;
//tacka koordinatnog pocetka za prvo i drugo stanje koordinatnog sistema
static double x_1 = -1, y_1 = -2, z_1 = 6;
static double x_2 = 4, y_2 = 1, z_2 = -4;
//parametri koji zavise od vremena
static double alfa_t = 0, beta_t = 0, gamma_t = 0;
static double x_t = 0, y_t = 0, z_t = 0;
static Vector4d q_t;
//kvaternioni koji predstavljaju pocetnu i krajnju rotaciju sistema
static Vector4d q_1;
static Vector4d q_2;
//parametri protoka vremena
static double tm = 10;
static double t = 0;
static int timerActive = 0;
//parametar koji odredjuje koji se algoritam koristi
static int odgovor;
//matrica rotacije
static Matrix3d A;
double fi;

//callback funkcije
static void on_timer(int value);
static void on_keyboard(unsigned char key, int x, int y);
static void on_reshape(int width,int height);
static void on_display();

void Euler2A(double alfa, double beta, double gamma, Matrix3d &E2A){

	//rotacija oko z-ose nekog tela se zove skretanje-yaw(ψ), oko y-ose propinjanje-pitch (θ) i oko x valjanje-roll (φ)
	//kada zelimo da rotiramo za ojlerove uglove dobijemo formulu A = RX2(φ)RY1(θ)RZ(ψ) - kada gledamo kao sopstveni sistem(y1 i x2 su ose nekih novih koordinatnih sistema)
	//kada to malo sredimo dobijemo A = RZ(ψ)RY(θ)RX(φ) ovde x,y i z predstavljaju ose svetskog koordinatnog sistema

	double alfasin = sin(alfa);
	double alfacos = cos(alfa);

	double betasin = sin(beta);
	double betacos = cos(beta);

	double gammasin = sin(gamma);
	double gammacos = cos(gamma);

	Matrix3d RXalfa;
	Matrix3d RYbeta;
	Matrix3d RZgamma;

	RXalfa << 1,       0,        0, 
	          0, alfacos, -alfasin,
	          0, alfasin,  alfacos;

	RYbeta << betacos, 0, betasin,
	                0, 1,       0,
	         -betasin, 0, betacos;

	RZgamma << gammacos, -gammasin, 0,
	           gammasin,  gammacos, 0,
	                  0,         0, 1;


	E2A = RZgamma * RYbeta * RXalfa;

	return;
}
//funkcija 2. : AxisAngle[A] - vraca jedinicni vektor p = (px, py, pz) i ugao alfa(0, pi) tako da A = Rp(alfa)
void AxisAngle(Matrix3d &A, Vector3d &p, double &angle){

	//vrsimo provere da li je A*A.transpose() = E i da li je A.det = 1
	Matrix3d E;

	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	Matrix3d E_test;

	E_test = A.transpose() * A;

	Matrix3d E_test_2;
	E_test_2 << abs(round(E_test(0, 0))), abs(round(E_test(0, 1))), abs(round(E_test(0, 2))),
	            abs(round(E_test(1, 0))), abs(round(E_test(1, 1))), abs(round(E_test(1, 2))),
	            abs(round(E_test(2, 0))), abs(round(E_test(2, 1))), abs(round(E_test(2, 2)));

	if(E_test_2 != E){
		std::cout << "Nazalost uneli ste matricu koja ne zadovoljava nase rigorozne kriterijume" << std::endl;
		return;
	}

	if(round(A.determinant()) != 1){
		std::cout << "Nazalost uneli ste matricu koja ne zadovoljava nase rigorozne kriterijume" << std::endl;
		return;
	}

	//racunamo vektor p - to je vektor oko koga se rotiramo tj. vektor koji se slika u samog sebe pri A
	//I nacin - znamo da jedna sopstvena vrednost matrice A mora biti 1 -> imamo formulu (A - 1*E)*p=0
	//odatle dobijamo 3 jednacine koje kad resimo dobijemo vektor p
	//II nacin - koji cemo i da radimo jeste da primetimo da su dva vektora iz (A - 1*E) normalna na p
	// iz toga sledi da iz njihovog vektorskog proizvoda mozemo dobiti p (pri tome ne biramo lin. zav. vektore) 
    Vector3d v1;

	v1 << A(0, 0)-1, A(0, 1), A(0, 2);

	Vector3d v2;

	v2 << A(1, 0), A(1, 1)-1, A(1, 2);

    //nomralizujemo ih
	v1 = v1.normalized();

	v2 = v2.normalized();

	//hvala andjo
	Vector3d zeroVector;

	zeroVector << 0, 0, 0;

	Vector3d crossy;
	crossy = v1.cross(v2);

	crossy(0) = round(crossy(0));
	crossy(1) = round(crossy(1));
	crossy(2) = round(crossy(2));

	//provera da li su lin zavisni vektori
	//podelimo svaku koordinato sa svakom - ako se dobije isti broj linearno su zavisni	
	if(crossy == zeroVector){
		v2 << A(2, 0), A(2, 1), A(2, 2)-1;
		v2 = v2.normalized();
	}


	p = v1.cross(v2);

	//normalizujemo i p
	p = p.normalized();

	//sad kad imamo p treba jos samo da nadjemo ugao rotacije
	//imamo vektor p, nadjemo neki vektor normalan na njega(npr neki koji smo vec koristili v1 ili v2)
	//na taj vektor primenimo matricu A koja rotira vektor v1 za ugao koji trazimo
	//sad samo primenimo formulu koja vraca ugao izmedju dva vektora i nas zadatak je gotov

	Vector3d v1p;

	v1p = A * v1;

	angle = std::acos(v1.dot(v1p));

	//na kraju samo proverimo koristeci mesoviti proizvod da li smo izabrali korektan smer prave p
	Matrix3d mesovita;

	mesovita << p, v1, v1p;

	if(mesovita.determinant() < 0){
		p = -1 * p;
	}

	return;
}
//funkcija 3 : Rodrigez[p, φ] - ulaz: prava p oko koje se rotira i ugao φ za koji se rotira; izlaz: matrica rotacije oko orjentisanog (jedinicnog) vektora p za ugao φ.
Matrix3d Rodrigez(Vector3d &p, double angle){

    //matrica koju vracamo se dobija iz ove formule:
	//Rp(φ) = ppT + cos φ (E − ppT) + sin φ p×

	//gde je p× matrica vektorskog mnozenja jedinicnim vektorom prave p

	Matrix3d Px;

	Px <<        0, -p(2, 0),  p(1, 0),
	       p(2, 0),        0, -p(0, 0),
	      -p(1, 0),  p(0, 0), 0;

	//jedinicna matrica
	Matrix3d E;

	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	//trazena matrica
	Matrix3d RodMat;

	RodMat = p*p.transpose() + cos(angle)*(E - p*p.transpose()) + sin(angle)*Px;

	return RodMat;
}
//Funkcija 4 : A2Euler[A] - za datu ortogonalnu matricu A, det A = 1, vraca Ojlerove uglove φ, θ, ψ
void A2Euler(Matrix3d &RodMat){


	double psi, teta, fi; 
	//Ulaz: Ortogonalna matrica A = (aij), det A = 1.
	//Izlaz: Ojlerovi uglovi ψ, θ, φ, A = Rz (ψ)Ry (θ)Rx (φ)

		if (RodMat(2, 0) < 1){
			if (RodMat(2, 0) > -1){ // jedinstveno resenje
				psi = atan2(RodMat(1, 0), RodMat(0, 0));
				teta = asin(-RodMat(2, 0));
				fi = atan2(RodMat(2, 1), RodMat(2, 2));
			} else { // nije jedinstveno, slucaj Ox3 = −Oz
				psi = atan2(-RodMat(0, 1), RodMat(1, 1));
				teta = M_PI/2;
				fi = 0;
			}
			} else { // nije jedinstveno, slucaj Ox3 = Oz
				psi = atan2(-RodMat(0, 1), RodMat(1, 1));
				teta = -M_PI/2;
			    fi = 0;
			}

	return;
}
//5. AxisAngle2Q[p, φ] - vraca jednicni kvaternion q = (x, y, z, w) tako da
//Cq = Rp(φ). Vektor p je jednicni.
void AxisAngle2Q(Vector3d& p, double angle, Vector4d& q){

	q << p(0, 0)*sin(angle/2), p(1, 0)*sin(angle/2), p(2, 0)*sin(angle/2), cos(angle/2);

	return;
}
//6. Q2AxisAngle[q] - vraca jedinicni vektor p = (px, py, pz) i ugao φ pripada [0, π]
//tako da kvaternion predstavlja rotaciju Rp(φ), tj. Cq = Rp(φ).
double QtoAxisAngle(Vector4d& q, Vector3d& vekP){

	double angleFiTmp = acos(q(3, 0));

	vekP << q(0, 0)/sin(angleFiTmp/2), q(1, 0)/sin(angleFiTmp/2), q(2, 0)/sin(angleFiTmp/2);

	angleFiTmp = 2*angleFiTmp;

	vekP = vekP.normalized();

	return angleFiTmp;
}


//cela poenta ovog domaceg
void SLERP(){
    
    //LERP in SLERP
    if(fi < M_PI/12 && fi > -M_PI/12){
        q_t = (1 - t/tm)*q_1 + (t/tm)*q_2;
        q_t = q_t.normalized();
    }
    else{
    q_t = sin(fi*(1 - t/tm))/sin(fi)*q_1 + sin(fi*(t/tm))/sin(fi)*q_2;
    }

    x_t = (1 - t/tm)*x_1 + (t/tm)*x_2;
    y_t = (1 - t/tm)*y_1 + (t/tm)*y_2;
    z_t = (1 - t/tm)*z_1 + (t/tm)*z_2;

    return;   
}
//Koordinatni sistem
void nacrtaj_KS(){

    GLUquadricObj *quadObj = gluNewQuadric();
	
    glPushMatrix();
    glColor3f(0, 1, 0.5);
    glTranslatef(3, 0, 0);
    glRotatef(90, 0, 1, 0);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0.25, 0.25);
    glTranslatef(0, 3, 0);
    glRotatef(-90, 1, 0, 0);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0, 1);
    glTranslatef(0, 0, 3);

    gluCylinder(quadObj, 0.1, 0, 0.3, 30, 30);
    glPopMatrix();

    glColor3f(0.25, 0.25, 1);
    glutSolidSphere(0.1, 30, 30);

    glBegin(GL_LINES);
        glColor3f(0, 1, 0.5);
        glVertex3d(0, 0, 0);
        glVertex3d(3, 0, 0);
        
        glColor3f(1, 0.25, 0.25);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 3, 0);
        
        glColor3f(1, 0, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 0, 3);
    glEnd();
}

void nacrtaj_pocetak_i_kraj(){
    Matrix3d E2A_1;

    Euler2A(alfa_1, beta_1, gamma_1, E2A_1);

    GLdouble matrix_1[16] = {E2A_1(0, 0), E2A_1(1, 0), E2A_1(2, 0), 0,
                          E2A_1(0, 1), E2A_1(1, 1), E2A_1(2, 1), 0,
                          E2A_1(0, 2), E2A_1(1, 2), E2A_1(2, 2), 0,
                          x_1, y_1, z_1, 1 };

    glPushMatrix();

    glMultMatrixd(matrix_1);
    nacrtaj_KS();

    glPopMatrix();

    Matrix3d E2A_2;

    Euler2A(alfa_2, beta_2, gamma_2, E2A_2);

    GLdouble matrix_2[16] = {E2A_2(0, 0), E2A_2(1, 0), E2A_2(2, 0), 0,
                          E2A_2(0, 1), E2A_2(1, 1), E2A_2(2, 1), 0,
                          E2A_2(0, 2), E2A_2(1, 2), E2A_2(2, 2), 0,
                          x_2, y_2, z_2, 1 };

    glPushMatrix();

    glMultMatrixd(matrix_2);
    nacrtaj_KS();

    glPopMatrix();

}

void izracunaj_kvaternione(){
    //calculating q1
    Matrix3d E2A_1;

    Euler2A(alfa_1, beta_1, gamma_1, E2A_1);

    Vector3d p_1;
    double angle_1;

    AxisAngle(E2A_1, p_1, angle_1);

    AxisAngle2Q(p_1, angle_1, q_1);

    //calculating q2
    Matrix3d E2A_2;

    Euler2A(alfa_2, beta_2, gamma_2, E2A_2);

    Vector3d p_2;
    double angle_2;

    AxisAngle(E2A_2, p_2, angle_2);

    AxisAngle2Q(p_2, angle_2, q_2);

    return;
}

int main(int argc, char* argv[]){

    std::cout << "unesite 1 ako zelite SLERP algoritam" << std::endl;
    std::cout << "unesite 2 ako (iz nekog razloga) zelite DERP algoritam" << std::endl;
    std::cin >> odgovor;
  
    //ucitavanje pozicija koordinatnih sistema
    //treba nam 3 ojlerova ugla i tacka centra novog koordinatnog sistema
    std::cout << "unesite alfa beta i gamma za prvo stanje Koordinatnog sistema:" << std::endl;
    std::cin >> alfa_1 >> beta_1 >> gamma_1;
    std::cout << "a: " << alfa_1 << " b: " << beta_1 << " g: " << gamma_1 << std::endl;
    alfa_1 = alfa_1*M_PI/180;
    beta_1 = beta_1*M_PI/180;
    gamma_1 = gamma_1*M_PI/180;

    //unos centra koordinatnog sistema
    std::cout << "unesite centar za prvo stanje koordiantnog sistema: " << std::endl;
    std::cin >> x_1 >> y_1 >> z_1;
    std::cout << "x: " << x_1 << " y: " << y_1 << " z: " << z_1 << std::endl;


    

    //sve isto samo za zavrsno stanje koordinatnog sistema
    std::cout << "unesite alfa beta i gamma za zavrsno stanje Koordinatnog sistema:" << std::endl;
    std::cin >> alfa_2 >> beta_2 >> gamma_2;
    std::cout << "a: " << alfa_2 << " b: " << beta_2 << " g: " << gamma_2 << std::endl;
    alfa_2 = alfa_2*M_PI/180;
    beta_2 = beta_2*M_PI/180;
    gamma_2 = gamma_2*M_PI/180;

    std::cout << "unesite centar za prvo stanje koordiantnog sistema: " << std::endl;
    std::cin >> x_2 >> y_2 >> z_2;
    std::cout << "x: " << x_2 << " y: " << y_2 << " z: " << z_2 << std::endl;

    std::cout << "unesite duzinu trajanja animacije: " << std::endl;
    std::cin >> tm;
    std::cout << "tm: " << tm << std::endl;

    izracunaj_kvaternione();

    double fi = acos(q_1(0)*q_2(0) + q_1(1)*q_2(1) + q_1(2)*q_2(2) + q_1(3)*q_2(3));
    if(fi > M_PI/2 || fi < -M_PI/2){
        q_1 = -q_1;
        fi = acos(q_1(0)*q_2(0) + q_1(1)*q_2(1) + q_1(2)*q_2(2) + q_1(3)*q_2(3)); 
    }
 
    //inicijalizuje se glut
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

    //Kreira se prozor
    glutInitWindowPosition(window_x, window_y);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("SLERPYDERP");

    //registruju se callback funkcije
    glutKeyboardFunc(on_keyboard);
    glutReshapeFunc(on_reshape);
    glutDisplayFunc(on_display);

    //postavljamo boju za ciscenje i omogucavamo test dubine
    glClearColor(0.87, 0.87, 0.87, 0);
    glEnable(GL_DEPTH_TEST);

    //podebljao sam linije
    glLineWidth(3);

    //program ulazi u glavnu petlju
    glutMainLoop();

    return 0;
}

static void on_keyboard(unsigned char key, int x, int y){
    
    switch(key){
        case 27:
            exit(0);
            break;
        case 'g':
        case 'G':
            // Pokrece se simulacija
            if (!timerActive) {
                glutTimerFunc(50, on_timer, 0);
                timerActive = 1;
            }
            break;
        case 's':
        case 'S':
            
            // Zaustavlja se simulacija
            timerActive = 0;
            break;
        case 'r':
        case 'R':            
            // Zaustavlja se simulacija
            timerActive = 0;
            t = 0;
            break;
    }
}


static void on_timer(int value)
{
    // Proverava se da li callback dolazi od odgovarajuceg tajmera
    if (value != 0)
        return;

    // Azurira se vreme simulacije
    t += 0.1;


    //proveravamo da li je vreme da zavrsimo animaciju
    if(t >= tm){
        t = 0;
        timerActive = 0;
        glutPostRedisplay();
        return;
    }

    // Forsira se ponovno iscrtavanje prozora
    glutPostRedisplay();

    // Po potrebi se ponovo postavlja tajmer
    if (timerActive)
        glutTimerFunc(50, on_timer, 0);
}
//menjanje dimenzija prozora
static void on_reshape(int width,int height){
    window_height = height;
    window_width = width;
}
//prikaz animacije
static void on_display(){

    // Brise se prethodni sadrzaj buffera
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Podesava se viewport
    glViewport(0, 0, window_width, window_height);

    // Podesava se projekcija
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, (float) window_width / window_height, 1, 200);

    // Podesava se pozicija sinteticke kamere
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(10, 10, 5,
             0, 0, 0,
             0, 1, 0);

    //ako nam zatreba provera kolika nam je matrica
    //GLfloat matrix[16]; 
    //glGetFloatv (GL_MODELVIEW_MATRIX , matrix);

    //pocetni i krajnji Koordinatni sistem
    nacrtaj_pocetak_i_kraj();

    glPushMatrix();

    if(odgovor == 1){

    //dobijamo trenutnu rotaciju u formi kvaterniona
    SLERP();
    
    //pretvaramo u dva koraka kvaternion u matricu rotacije
    Vector3d vekP;
    double angle = QtoAxisAngle(q_t, vekP);

    A = Rodrigez(vekP, angle);

    }
    else if(odgovor == 2){

    alfa_t = (1 - t/tm)*alfa_1 + t/tm*alfa_2;
    beta_t = (1 - t/tm)*beta_1 + t/tm*beta_2;
    gamma_t = (1 - t/tm)*gamma_1 + t/tm*gamma_2;
    x_t = (1 - t/tm)*x_1 + t/tm*x_2;
    y_t = (1 - t/tm)*y_1 + t/tm*y_2;
    z_t = (1 - t/tm)*z_1 + t/tm*z_2;

    Euler2A(alfa_t, beta_t, gamma_t, A);

    }
    else{
        std::cout << "ako ne mozes da pratis jednostavne instrukcije verovatno ne bi ni skapirao sta se desava ovde..." << std::endl;
        exit(1);
    }

    //matrica kojom predstavljamo rotacije i translacije
    GLdouble matrix[16] = {
                           A(0, 0), A(1, 0), A(2, 0), 0,
                           A(0, 1), A(1, 1), A(2, 1), 0,
                           A(0, 2), A(1, 2), A(2, 2), 0,
                               x_t,     y_t,     z_t, 1 
                          };

    glMultMatrixd(matrix);

    nacrtaj_KS();

    glPopMatrix();

    //Svetski Koordinatni sistem 
    glColor3f(0.25, 0.25, 1);
    glutSolidSphere(0.1, 30, 30);

    glBegin(GL_LINES);
        glColor3f(0, 1, 0.5);
        glVertex3d(-200, 0, 0);
        glVertex3d(200, 0, 0);
        
        glColor3f(1, 0.25, 0.25);
        glVertex3d(0, -200, 0);
        glVertex3d(0, 200, 0);
        
        glColor3f(1, 0, 1);
        glVertex3d(0, 0, -200);
        glVertex3d(0, 0, 200);
    glEnd();

    glutSwapBuffers();


}
