#include "alg.hpp"
#include <iostream>

using namespace std;

#define PI 3.14159265

vec2 mid;

template <typename T> void swap(vector<T>* vec, int indexA, int indexB) {
    if (indexA == indexB) return;
    T temp = vec->at(indexA);
    vec->at(indexA) = vec->at(indexB);
    vec->at(indexB) = temp;
}

int partition(vector<double>* vec, int startIndex, int endIndex) {
    int pivot = vec->at(endIndex);
    int i = startIndex - 1;
    for (int j = startIndex; j < endIndex; j++) {
        if (vec->at(j) < pivot) {
            i += 1;
            swap(vec, i, j);
        }
    }
    i += 1;
    swap(vec, i, endIndex);
    return i;
}


void quickSortStep(vector<double>* vec, int startIndex, int endIndex) {
    if (endIndex <= startIndex) return; // recursion f0 case

    int pivotIndex = partition(vec, startIndex, endIndex);
    quickSortStep(vec, startIndex, pivotIndex - 1);
    quickSortStep(vec, pivotIndex + 1, endIndex);
}


void quickSort(vector<double>* vec) {
    quickSortStep(vec, 0, vec->size() - 1);
}


void merge(vector<double>* vec, vector<double>* leftVec, vector<double>* rightVec) {
    vector<double> vecCopy;
    int size = vec->size();
    int leftSize = size / 2;
    int rightSize = size - leftSize;
    int l = 0;
    int r = 0;
    int step = 0;
    while (l < leftSize && r < rightSize) {
        if (leftVec->at(l) < rightVec->at(r)) {
            vec->at(step) = leftVec->at(l);
            l += 1;
        }
        else {
            vec->at(step) = rightVec->at(r);
            r += 1;
        }
        step++;
    }
    while (l < leftSize) {
        vec->at(step) = leftVec->at(l);
        l += 1;
        step++;
    }
    while (r < rightSize) {
        vec->at(step) = rightVec->at(r);
        r += 1;
        step++;
    }
}


void mergeSort(vector<double>* vec) {
    int subVecSize = vec->size();
    if (subVecSize < 2) return; //already sorted

    //two subvectors split in half
    int half = subVecSize / 2; // int division = floor

    vector<double> leftSubVec;
    for (int i = 0; i <= half - 1; i++) {
        leftSubVec.push_back(vec->at(i));
    }
    vector<double> rightSubVec;
    for (int i = half; i < subVecSize; i++) {
        rightSubVec.push_back(vec->at(i));
    }
    mergeSort(&leftSubVec);
    mergeSort(&rightSubVec);
    merge(vec, &leftSubVec, &rightSubVec);
}


int partitionVec2(vector<vec2>* points, int startIndex, int endIndex, const char& parameter) {
    if (parameter == 'x') {
        double pivot = points->at(endIndex).x;
        int i = startIndex - 1;
        for (int j = startIndex; j < endIndex; j++) {
            if (points->at(j).x < pivot) {
                i += 1;
                swap(points, i, j);
            }
        }
        i += 1;
        swap(points, i, endIndex);
        return i;
    }
    else if (parameter == 'y') {
        double pivot = points->at(endIndex).y;
        int i = startIndex - 1;
        for (int j = startIndex; j < endIndex; j++) {
            if (points->at(j).y < pivot) {
                i += 1;
                swap(points, i, j);
            }
        }
        i += 1;
        swap(points, i, endIndex);
        return i;
    }
    else {
        throw runtime_error("Invalid Parameter!");
    }
}

void quickSortStepVec2(vector<vec2>* vec, int startIndex, int endIndex, const char& parameter) {
    if ((parameter == 'x') or (parameter == 'y')) {
        if (endIndex <= startIndex) return; // recursion f0 case

        int pivotIndex = partitionVec2(vec, startIndex, endIndex, parameter);
        quickSortStepVec2(vec, startIndex, pivotIndex - 1, parameter);
        quickSortStepVec2(vec, pivotIndex + 1, endIndex, parameter);
    }
    else {
        throw runtime_error("Invalid Parameter!");
    }
}

void quickSortVec2x(vector<vec2>* points) {
    quickSortStepVec2(points, 0, points->size() - 1, 'x');
}

void quickSortVec2y(vector<vec2>* points) {
    quickSortStepVec2(points, 0, points->size() - 1, 'y');
}


double orientation(const vec2& p1, const vec2& p2, const vec2& p3) {
    //Compare the slope between two known points. So, in the end of the day, we can know the orientation.
    double res = (p3.y - p2.y) * (p2.x - p1.x) - (p2.y - p1.y) * (p3.x - p2.x);
    if (res > 0) {
        //CCW
        return 1;
    }
    else if (res < 0) {
        //CW
        return -1;
    }
    else {
        //Colinear
        return 0;
    }
}

double dist(const vec2& p1, const vec2& p2) {
    return sqrt(pow(p2.y - p1.y, 2) + pow(p2.x - p1.x, 2));
}

vector<vec2> jarvis(vector<vec2>* points) {

    vector<vec2> res{};

    //First Step -> Find P0
    quickSortVec2x(points);

    /*Pathological Case : If there are extreme points with the same X value,
    we will take the point with the highest Y value among these */
    int on_hull = 0;
    if (points->at(0).x == points->at(1).x) {
        for (int i = 1; i < (points->size() - 1); i++) {
            if (points->at(i).x != points->at(i + 1).x) {
                //Index 'i' is the last position which contains a candidate point.
                double y_max = points->at(0).y;
                for (int j = 1; j <= i; j++) {
                    if (points->at(j).y > y_max) {
                        y_max = points->at(j).y;
                        on_hull = j;
                    }
                }
                break;
            }
        }
    }

    //Now, we have p0!


    //Now, let's go to the magic!
    bool cond = true;

    while (cond) {
        res.push_back(points->at(on_hull));
        vec2 next_point = points->at(0); //The fist point in the set
        vec2 point_in_the_hull = points->at(on_hull);
        for (int i = 0; i < points->size(); i++) {
            /* 1- point_in_the_hull is the last point appended in the convex hull, and we are sure
            that this point belongs to the convex hull subset.
               2- next_point is the next candidate to be appended. We must check the orientation!
               3- points->at(i) is used to run the set of points */
            double o = orientation(point_in_the_hull, next_point, points->at(i));
            /*We update the next_point when the orientation is CCW(o == 1).But also, when the actual
            next_point is the last one we added to the convex hull, and sure, we shouldn't add a point
            which is already in the result subset.
            And, sure, if we have 3 collinear points, we must take the far way point to minimize
            the result subset-> The point in that moment of the loop is the new candidate*/
            if ((next_point == point_in_the_hull) or (o == 1) or (o == 0 and
                dist(point_in_the_hull, points->at(i)) >= dist(point_in_the_hull, next_point))) {
                next_point = points->at(i);
                on_hull = i;
            }
        }
        //Now, we have the right candidate.

        //Loop Conditional
        if (points->at(on_hull) == res.at(0)) {
            cond = false;
        }
    }

    res.push_back(res.at(0));
    return res;
}

//MergeHull

/*To determine the quadrant
(+, +) -> Fist Quadrant
(-, +) -> Second Quadrant
(-, -) -> Third Quadrant
(+, -) -> Forth Quadrant */
int quad(vec2 p) {
    if (p.x >= 0 && p.y >= 0)
        return 1;
    if (p.x <= 0 && p.y >= 0)
        return 2;
    if (p.x <= 0 && p.y <= 0)
        return 3;
    return 4;
}

//Used inside std::sort function -> Sorting in CCW order
bool compare(vec2 p1, vec2 q1) {
    vec2 p(p1.x - mid.x, p1.y - mid.y);
    vec2 q(q1.x - mid.x, q1.y - mid.y);

    int one = quad(p);
    int two = quad(q);

    if (one != two)
        return (one > two);
    return (p.y * q.x > q.y * p.x);
}

// Checks the orientation of the triplet (a, b, c)
// 0 -> a, b and c are collinear
// 1 -> Clockwise
// -1 -> Counterclockwise
int check_ori(vec2 a, vec2 b, vec2 c) {
    double res = (b.y - a.y) * (c.x - b.x) - (c.y - b.y) * (b.x - a.x);

    if (res == 0)
        return 0;
    if (res > 0)
        return 1;
    return -1;
}

vector<vec2> subproblem_jar(vector<vec2>* points) {
    vector<vec2> res{};

    //First Step -> Find P0
    quickSortVec2x(points);

    /*Pathological Case : If there are extreme points with the same X value,
    we will take the point with the highest Y value among these */
    int on_hull = 0;
    if (points->at(0).x == points->at(1).x) {
        for (int i = 1; i < (points->size() - 1); i++) {
            if (points->at(i).x != points->at(i + 1).x) {
                //Index 'i' is the last position which contains a candidate point.
                double y_max = points->at(0).y;
                for (int j = 1; j <= i; j++) {
                    if (points->at(j).y > y_max) {
                        y_max = points->at(j).y;
                        on_hull = j;
                    }
                }
                break;
            }
        }
    }

    //Now, we have p0!


    //Now, let's go to the magic!
    bool cond = true;

    while (cond) {
        res.push_back(points->at(on_hull));
        vec2 next_point = points->at(0); //The fist point in the set
        vec2 point_in_the_hull = points->at(on_hull);
        for (int i = 0; i < points->size(); i++) {
            /* 1- point_in_the_hull is the last point appended in the convex hull, and we are sure
            that this point belongs to the convex hull subset.
               2- next_point is the next candidate to be appended. We must check the orientation!
               3- points->at(i) is used to run the set of points */
            double o = orientation(point_in_the_hull, next_point, points->at(i));
            /*We update the next_point when the orientation is CCW(o == 1).But also, when the actual
            next_point is the last one we added to the convex hull, and sure, we shouldn't add a point
            which is already in the result subset.
            And, sure, if we have 3 collinear points, we must take the far way point to minimize
            the result subset-> The point in that moment of the loop is the new candidate*/
            if ((next_point == point_in_the_hull) or (o == 1) or (o == 0 and
                dist(point_in_the_hull, points->at(i)) >= dist(point_in_the_hull, next_point))) {
                next_point = points->at(i);
                on_hull = i;
            }
        }
        //Now, we have the right candidate.

        //Loop Conditional
        if (points->at(on_hull) == res.at(0)) {
            cond = false;
        }
    }

    // Sorting the points in the anti-clockwise order
    mid = { 0, 0 };
    int n = res.size();
    for (int i = 0; i < n; i++) {
        mid.x += res[i].x;
        mid.y += res[i].y;
    }
    mid.x = mid.x / n;
    mid.y = mid.y / n;
    //Now we have the centroid!
    sort(res.begin(), res.end(), compare);

    return res;
}

vector<vec2> mergeHullMerger(vector<vec2> a, vector<vec2> b) {
    // n1 -> number of points in polygon a
    // n2 -> number of points in polygon b
    int n1 = a.size(), n2 = b.size();

    int ia = 0, ib = 0;
    for (int i = 1; i < n1; i++)
        if (a[i].x > a[ia].x)
            ia = i;

    // ib -> leftmost point of b
    for (int i = 1; i < n2; i++)
        if (b[i].x < b[ib].x)
            ib = i;

    // finding the upper tangent
    int inda = ia, indb = ib;
    bool done = 0;
    while (!done) {
        done = 1;
        while (check_ori(b[indb], a[inda], a[(inda + 1) % n1]) > 0)
            inda = (inda + 1) % n1;

        while (check_ori(a[inda], b[indb], b[(n2 + indb - 1) % n2]) < 0) {
            indb = (n2 + indb - 1) % n2;
            done = 0;
        }
    }

    int uppera = inda, upperb = indb;

    inda = ia, indb = ib;
    done = 0;
    while (!done) { // finding the lower tangent
        done = 1;
        while (check_ori(a[inda], b[indb], b[(indb + 1) % n2]) > 0)
            indb = (indb + 1) % n2;

        while (check_ori(b[indb], a[inda], a[(n1 + inda - 1) % n1]) < 0) {
            inda = (n1 + inda - 1) % n1;
            done = 0;
        }
    }

    int lowera = inda, lowerb = indb;
    
    vector<vec2> ret;

    // ret contains the convex hull after merging the two convex hulls
    // with the points sorted in anti-clockwise order
    int ind = uppera;
    ret.push_back(a[uppera]);
    while (ind != lowera) {
        ind = (ind + 1) % n1;
        ret.push_back(a[ind]);
    }

    ind = lowerb;
    ret.push_back(b[lowerb]);
    while (ind != upperb) {
        ind = (ind + 1) % n2;
        ret.push_back(b[ind]);
    }
    return ret;
}

// Returns the convex hull for the given set of points
vector<vec2> mergeHullDivide(vector<vec2> a) {
    // If the number of points is less than 6 then the
    // function uses the jarvis algorithm to find the
    // convex hull
    if (a.size() <= 5) {
        return subproblem_jar(&a);
    }

    // left contains the left half points
    // right contains the right half points
    vector<vec2> left, right;
    for (int i = 0; i < a.size() / 2; i++)
        left.push_back(a[i]);
    for (int i = a.size() / 2; i < a.size(); i++)
        right.push_back(a[i]);

    // convex hull for the left and right sets
    vector<vec2> left_hull = mergeHullDivide(left);
    vector<vec2> right_hull = mergeHullDivide(right);

    // merging the convex hulls
    return mergeHullMerger(left_hull, right_hull);
}

vector<vec2> mergeHull(vector<vec2> points) {

    //First, we need to sort according to the x-axis
    quickSortVec2x(&points);

    //Do the Magic!
    vector<vec2> res = mergeHullDivide(points);
    res.push_back(res.at(0));

    return res;
}

//Advancing Front

vector<vec2> advancingFront(vector<vec2>& inputPoints) {
    
    vector<vec2> convexHull = jarvis(&inputPoints);
    double x = 0, y = 0;
    for (auto &point : convexHull) {
        x = x + point.x;
        y = y + point.y;
    }
    // - 1 --> Because point[0] appears twice.
    x = x / (convexHull.size() - 1);
    y = y / (convexHull.size() - 1);

    mid.x = x; mid.y = y;
    //sort(convexHull.begin(), convexHull.end(), compare);

    vector<pair<vec2, vec2>> boundary{};
    for (size_t i = 0; i < convexHull.size() - 1; i++) {
        vec2 a = convexHull.at(i);
        vec2 b = convexHull.at(i + 1);
        boundary.emplace_back(std::make_pair(a, b));
    }

    return adf_magic(inputPoints, boundary);
}

bool intersec2D(const pair<vec2, vec2>& r, const pair<vec2, vec2>& s) {

    vec2 k = r.first; vec2 l = r.second;
    vec2 m = s.first; vec2 n = s.second;
    
    if (k == l || m == n) {
        std::cout << "Problema 1 em intersec2D" << std::endl;
        return false;
    }
    
    if ((k == m && l == n) || (k == n && l == m)) {
        std::cout << "Segmentos iguais foram passados para intersec2D\n";
        return false;
    }

    double det = (n.x - m.x) * (l.y - k.y) - (n.y - m.y) * (l.x - k.x);

    if (det == 0.0) { return false; } // não há intersecção
    else { return true; } //há intersecção

}

static bool rotationIndexPosition(const vector<pair<vec2, vec2>>* points, const vec2& q) {
    //Points is a vector of polygon's vertices, e.g, for a square, we have p0, p1, p2, p3.

    vector<double> angles{};

    for (int i = 0; i < points->size(); i++) {
        //We'll need -> vec2(p_i - q) and (P_i+1 - q) -> We gonna call vec2(a) and vec2(b)
        //Extra Case: In the last loop, we gonna have -> vec(p_max - q) and vec2(p_0 - q)

        vec2 a{};
        vec2 b{};

        a = points->at(i).first - q;
        b = points->at(i).second - q;

        //The angle between a and b = arc cos((a . b) / (|a| |b|))
        double teta = acos((a.dot(b)) / (a.mag() * b.mag()));
        teta = teta * 180.0 / PI; //Rad -> Degrees

        angles.push_back(teta);
    }

    double sum = 0;

    for (int i = 0; i < angles.size(); i++) {
        sum += angles.at(i);
    }

    //I will tolerate an error of up to 1e-3
    //+-1 -> Inside -> True
    if (((1 <= sum / 360) and (sum / 360 <= 1 + 1e-3)) or ((-1 >= sum / 360) and (sum / 360 >= -1 - 1e-3))) {
        return true;
    }
    else {
        return false;
    }
}

static bool isConvexHull(const vector<pair<vec2, vec2>>& originalBoundary, const vec2& point) {
    for (auto &edge: originalBoundary) {
        if ((edge.first == point) || (edge.second == point)) {
            return true;
        }
    }
    return false;
}

vector<vec2> adf_magic(vector<vec2> inputPoints, vector<pair<vec2, vec2>>& boundary){
    
    vector<pair<vec2, vec2>> originalBoundary{ boundary.begin(), boundary.end() };

    //Empty Final Set
    vector<vec2> resultado{};

    //Let's put the boundary in a stack
   
    cout << "Convex Hull:" << endl;
    for (auto& edge: boundary) {
        std::cout << "(" << edge.first.x << ", " << edge.first.y << "), (" <<
            edge.second.x << ", " << edge.second.y << "),\n";
    }
    
    std::cout << "Tamanho inicial da Pilha: " << boundary.size() << std::endl;

    bool stop_condition = false; int main_while_control = 1;
    vector<pair<vec2, double>> rank_distance{};
    //TOPO DA PILHA = PRIMEIRO ELEMENTO DO VECTOR
    while (!stop_condition) {

        std::cout << "LOOP PRINCIPAL: " << main_while_control << endl;

        pair<vec2, vec2> current_edge = boundary.front();
        vec2 a = current_edge.first; vec2 b = current_edge.second;
        
        if (a == b) {
            boundary.erase(boundary.begin());
            std::cout << "Pulei iteration pois a == b\n";
            continue;
        }

        for (auto &point: inputPoints) {
            if(!(point == a) && !(point == b)){
                double determinante = (point.x - a.x) * (b.y - a.y) - (point.y - a.y) * (b.x - a.x);
                if (determinante > 0.0) {
                    if (rotationIndexPosition(&boundary, point) || isConvexHull(boundary, point)) {
                        double distancia_a = sqrt(std::pow(point.x - a.x, 2) + std::pow(point.y - a.y, 2));
                        double distancia_b = sqrt(std::pow(point.x - b.x, 2) + std::pow(point.y - b.y, 2));
                        double distancia = (distancia_a + distancia_b) / 2;
                        rank_distance.emplace_back(std::make_pair(point, distancia));
                    }
                }
            }
        }

        auto sorting_dist = [](const pair<vec2, double>& a, const pair<vec2, double>& b) {
            return a.second < b.second;
        };
        std::sort(rank_distance.begin(), rank_distance.end(), sorting_dist);

        std::cout << "Tamanho do rank distance: " << rank_distance.size() << std::endl;

        //TODO: Analise para evitar arestas que se cruzem

        vec2 added_point{};
        if (rank_distance.size() == 0 || resultado.size() / 3 == 150) {
            //std::cout << "RANK DISTANCE = 0!!!\n";
            std::cout << "Boundary final:\n";
            for (auto &edge : boundary) {
                std::cout << "(" << edge.first.x << ", " << edge.first.y << "), (" <<
                    edge.second.x << ", " << edge.second.y << "),\n";
            }
            return resultado;
        }
        added_point = rank_distance.front().first;
        std::cout << "Target Final: " << added_point << std::endl;

        resultado.push_back(a);
        resultado.push_back(added_point);
        resultado.push_back(b);
                
        std::cout << "Numero de Triangulos neste momento: " << resultado.size() / 3 << endl;

        //Limpando o rank_distance para a proxima iteração
        rank_distance.clear();

        //Now we need to advance the front
        boundary.erase(boundary.begin()); //100% sure b_a it was already at the front.
        
        //target_b
        int aux = 0;
        bool other = true;
        for (auto& edge : boundary) {
            if ((edge.first == added_point && edge.second == b) || (edge.first == b && edge.second == added_point)) {
                boundary.erase(boundary.begin() + aux);
                other = false;
                break;
            }
            aux = aux + 1;
        }

        if (other) {
            std::cout << "target_b adicionado a boundary!\n";
            boundary.insert(boundary.begin(), std::make_pair(added_point, b));
        }

        //a_target or target_a are at boundary?
        aux = 0;
        bool boolean = true;
        for (auto& edge : boundary) {
            if ((edge.first == a && edge.second == added_point) || (edge.first == added_point && edge.second == a)) {
                boundary.erase(boundary.begin() + aux);
                boolean = false;
                break;
            }
            aux = aux + 1;
        }

        if (boolean) {
            std::cout << "a_target adicionado a boundary!\n";
            boundary.insert(boundary.begin(), std::make_pair(a, added_point));
        }

        main_while_control++;

        if ((boundary.size() > 0) && !(boundary.back().second == boundary.front().first)) {
            boundary.back().second = boundary.front().first;
        }

        cout << "Boundary Atual:\n";
        for (auto &edge: boundary) {
            std::cout << "(" << edge.first.x << ", " << edge.first.y << "), (" <<
                edge.second.x << ", " << edge.second.y << "),\n";
        }

        //Stop condition
        if (boundary.empty()) {
            std::cout << "STOP CONDITION!!\n";
            stop_condition = true;
        }
    }

    return resultado;
}