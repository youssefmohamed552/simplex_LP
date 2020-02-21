#include <iostream>
// #include <bits/stdc++.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>


class Matrix{
  private:
    int m_m;
    int m_n;
    std::vector<std::vector<double>> m_matrix;
  public:
    Matrix(const int m, const int n)
      : m_m(m), m_n(n), m_matrix(m, std::vector<double>(n, 0)){}
    ~Matrix(){}
    int m() const { return m_m; }
    int n() const { return m_n; }

    void raw_reduce(Matrix& m,  int i, int j){
      // Matrix& m = *this;
      double num = m[i][j];
      for(int k = 0; k < this->n(); k++){
        if(m[i][k] != 0) m[i][k] /= num;
      }
      for(int k = i+1; k < this-> m(); k++){
        num = m[k][j];
        for(int l = 0; l < m.n(); l++){
          m[k][l] -= (m[i][l] * num);
        }
      }
      for(int k = i-1; k >= 0; k--){
        num = m[k][j];
        for(int l = 0; l < m.n(); l++){
          m[k][l] -= (m[i][l] * num);
        }
      }
    }

    Matrix& rref(){
      Matrix& rref = *this;
      for(int i = 0; i < rref.m(); i++){
        for(int j = 0; j < rref.n(); j++){
          if(rref[i][j] == 0) continue;
          // std::cout << rref;
          raw_reduce(rref,i,j);
          // std::cout << rref;
          break;
        }
      }
      return rref;
    }

    int rank(){
      Matrix temp = this->rref();
      // std::cout << "after rref\n";
      int rank = 0;
      for(int i = 0; i < temp.m(); i++){
        for(int j = 0; j < temp.n()-1; j++){
          if(temp[i][j] == 1){ 
            rank++; 
            break;
          }
        }
      }
      return rank;
    }

    bool is_full_rank(){
      Matrix temp = *(this);
      int r = temp.rank();
      // std::cout << "rank is : " << r << std::endl;
      return m_m == r;
    }

    std::vector<double> get_solution(){
      rank();
      // std::cout << *(this);
      std::vector<double> s(m_n - m_m - 1);
      for(int i = 0; i < m_m; i++){
        for(int j = 0; j < s.size(); j++){
          if(m_matrix[i][j] == 1){
            s[j] = m_matrix[i].back();
            break;
          }
        }
      }
      return s;
    }

    std::vector<double>& operator[](size_t i){ return m_matrix[i]; }
    friend std::ostream& operator<<(std::ostream& out,  Matrix& A){
      std::cout << "printing matrix\n";
      for(size_t i = 0; i < A.m(); i++){
        for(size_t j = 0; j < A.n(); j++){
          out << A[i][j] << " ";
        }
        out << std::endl;
      }
      return out;
    }
};

std::ostream& operator<<(std::ostream& out, const std::list<int>& l){
  for(auto e : l){
    out << e << " -> ";
  }
  out << "null" << std::endl;
  return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<double>& l){
  for(auto e : l){
    out << e << " -> ";
  }
  out << "null" << std::endl;
  return out;
}


class Simplex{
  private:
    Matrix m_input;
    Matrix m_augmented;
    std::vector<bool>m_sol;
  public:
    Simplex(Matrix& input)
      :m_input(input), m_augmented(Matrix(input.m(), input.n() + input.m())), m_sol(std::vector<bool>(input.m() + input.n() -1, true))
    {
      // int rank = input.rank();
      // if(input.is_full_rank()){
        // std::cout << "THE MATRIX IS FULL RANK\n";
      // }
      for(int i = 0; i < input.m(); i++){
        for(int j = 0; j < input.n() -1; j++){
          m_augmented[i][j] = input[i][j];
        }
      }
      for(int i = 0; i < input.m(); i++){
        // std::cout << "here\n";
        m_augmented[i][m_augmented.n()-1] = input[i][input.n()-1];
        m_augmented[i][m_input.n()+i-1] = 1;
      }
      // std::cout << m_augmented;
    }
    std::list<int> solution(){
      std::list<int> ans;
      for(int i = 0; i < m_sol.size(); i++){
        if(m_sol[i]) ans.push_back(i+1);
      }
      return ans;
    }

    bool is_feasible_solution(Matrix& m, std::vector<double>& x){
      // std::cout << x;
      // std::cout << m;
      for(int i = 0; i < m.m(); i++){
        double sum = 0.0;
        for(int j = 0; j < x.size(); j++){
          sum += (m[i][j] * x[j]);
          // std::cout << sum << " go ";
        }
        if(std::abs(sum - m[i].back()) > 0.0001) return false;
        // std::cout << sum << " = ";
      }
      // std::cout << "end \n";
      return true;
    }

    bool solve_basic_variables(std::list<int>& a){
      Matrix temp = m_augmented;
      for(auto j : a){
        for(int i = 0; i < m_augmented.m(); i++){
          temp[i][j] = 0;
        }
      }
      if(!temp.is_full_rank()){
        return false;
      }
      std::vector<double> sol = temp.get_solution();
      // std::cout << sol;
      // std::cout << m_input;
      if(!is_feasible_solution(m_input, sol)){
        // std::cout << "didn't pass\n";
        return false;
      }
      temp.rref();
      // std::cout << temp;
      bool feasible = true;
      for(int i = 0;i < temp.m(); i++){
        if(temp[i].back() < 0.0 ){
          // std::cout << "found an infeasible solution of value : " << temp[i].back() << " and at i = " << i << "\n";
          feasible = false;
          break;
        }
      }

      
      if(!feasible){
        // std::cout << "INFEASIBLE\n";
        return false;
      }
      // std::cout << temp;
      return true;
    }
    bool combination(std::vector<int>& t, std::list<int>& a,const int s, const int k){
      if(k == 0){
        // std::cout << a;
        bool t = solve_basic_variables(a); 
        if(t){
          // std::cout << "it's true for the combination " << a;
          for(auto e : a){
            m_sol[e] = false;
          }
        }
        return t;
      }
      int feasible = false;
      for(int i = s; i <= t.size() - k; i++){
        a.push_back(t[i]);
        feasible |= combination(t,a, i+1, k-1);
        a.pop_back();
      }
      return feasible;
    }

    bool solve(){
      if(!m_input.is_full_rank()){
        // std::cout << "INFEASIBLE\n";
        return false;
      }
      std::list<int> comp;
      std::vector<int> total(m_augmented.n()-1);
      for(int i = 0; i < total.size(); i++){ total[i] = i; }

      return combination(total, comp,0, m_augmented.n()-m_augmented.m()-1);
      
    }

};



int main(void){
  std::cout << "Welcome to Simplex Algorithm !!\n";
  int k = 0;
  std::cin >> k;
  for(int x = 0; x < k; x++){
    int n = 0, m = 0;
    std::cin >> n >> m;
    Matrix A(m , (n+1));
    for(int i = 0; i < m; i++){
      for(int j = 0; j <= n; j++){
        std::cin >> A[i][j];
      }
    }
    // std::cout << A;
    // clock_t start, end;
    // start = clock();
    Simplex simplex(A);
    if(!simplex.solve()){
      std::cout << "INFEASIBLE\n";
    }
    else{
      // std::cout << "FEASIBLE\n";
      // std::cout << simplex.solution();
      for(auto e : simplex.solution()){
        std::cout << e << " ";
      }
      std::cout << std::endl;
    }
    // end = clock();
    // double time_taken = double(end - start)/ (double(CLOCKS_PER_SEC));
    // std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5);
    // std::cout << " sec " << std::endl;
    // break;

  }
  return 1;
}
