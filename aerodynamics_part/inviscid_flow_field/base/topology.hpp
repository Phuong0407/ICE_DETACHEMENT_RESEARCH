#ifndef topology_hpp
#define topology_hpp

namespace geometry{

template<unsigned int topodim, typename index_t>
class topology {
private:
    index_t id;
    
public:
    topology() = default;
    explicit topology(index_t id) : id(id) {}
    
public:
    index_t id() const { return this->id; }
    static constexpr unsigned int get_topodim() { return topodim; }
};
    
} // end namespace geometry
    
#endif