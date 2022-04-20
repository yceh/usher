#include <istream>
#include <deque>
#include <cstdio>
#include <vector>
static void print_content(std::deque<char>& fifo){
    for(auto c:fifo){
        fputc(c, stdout);
    }
    fputc('\n', stdout);
}
int main(int argc, char** argv){
    auto pos_fh=fopen(argv[1], "r");
    std::deque<int> positions;
    while (true) {
        int pos;
        int ret=fscanf(pos_fh, "%d\n",&pos);
        if (ret==1) {
            positions.push_back(pos);
        }
        if (ret==EOF) {
            break;
        }
    }
    std::deque<char> fifo;
    int pos=0;
    while (true) {
        auto read=fgetc(stdin);
        pos++;
        while (positions.front()+80<pos) {
            positions.pop_front();
            if (positions.empty()) {
                return 0;
            }
        }
        if (read==EOF) {
            print_content(fifo);
            break;
        }
        if (fifo.size()==150) {
            if (positions.front()+74<pos&&pos<positions.front()+76) {
                print_content(fifo);        
            }
            fifo.pop_front();
        }
        fifo.push_back(read);
    }
}