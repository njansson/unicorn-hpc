#include <dolfin.h>
#include <unicorn/init.h>
#include <unicorn/util.h>

using namespace dolfin;
using namespace dolfin::unicorn;

void pre(Mesh& mesh) {
}

void post(Mesh& mesh) {
}

void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time) {

}

int main(int argc, char *argv[]) {
  timeval s_time;
  gettimeofday(&s_time, NULL);
  Mesh mesh;  
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;  

  /*
   * Initialize unicorn, load mesh, parse parameters etc
   */
  unicorn_init(argc, argv, mesh, chkp, w_limit, iter);
  dolfin_set("output destination", "silent");

  /*
   * Start solver iteration
   */
  unicorn_solve(mesh, chkp, w_limit, s_time, iter, &pre, &post, &solve);

  dolfin_finalize();
  return 0;
}
