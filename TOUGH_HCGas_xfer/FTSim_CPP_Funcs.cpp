#include <iostream>
#include <unordered_map>

template<class T>
struct Hash_Key_Comp
{
    bool operator()(T x, T y) const
    {
      return (x[0] == y[0]
           && x[1] == y[1]
           && x[2] == y[2]
           && x[3] == y[3]
           && x[4] == y[4]);
    }
};

struct Char_Hash
{
    int operator()(char *id ) const
    {
        int seed = 1311311;
        int hash = 0;

        hash = (hash * seed) + (id[0]);
        hash = (hash * seed) + (id[1]);
        hash = (hash * seed) + (id[2]);
        hash = (hash * seed) + (id[3]);
        hash = (hash * seed) + (id[4]);

        return hash & (0x7FFFFFFF);
    }
};

extern "C"
{
  void check_edges(int num_elem, int num_conx, int *n1, int *n2, char **name1,
                   char **name2, char **elems)
  {
    std::unordered_map<char *, int, Char_Hash, Hash_Key_Comp<char *> > name_map;
    name_map.reserve(num_elem);

    for (int i = 0; i < num_elem; ++i)
      name_map[elems[i]] = i+1;

#ifdef USE_OMP
#pragma omp for
#endif
    for (int i = 0; i < num_conx; ++i)
    {
      auto name1_lookup = name_map.find(name1[i]);
      auto name2_lookup = name_map.find(name2[i]);

      if (name1_lookup != name_map.end())
        n1[i] = (*name1_lookup).second;
      else
        std::cout << "Reference to unknown element "
                  << name1[0] << name1[1] << name1[2] << name1[3] << name1[4]
                  << " at connection number " << i+1;

      if (name2_lookup != name_map.end())
        n2[i] = (*name2_lookup).second;
      else
        std::cout << "Reference to unknown element "
                  << name1[0] << name1[1] << name1[2] << name1[3] << name1[4]
                  << " at connection number " << i+1;
    }

    return;
  }
}
