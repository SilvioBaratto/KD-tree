template<typename T>
class point{
    public:
        point(T x, T y):
            coords(x, y) {}

        T x(){return coords.first;}
        T y(){return coords.second;}
        T index(int idx){
            if(idx == 0) return coords.first;
            return coords.second;
        }

    private:
        std::pair<T, T> coords;
};
