#ifndef Materials_h
#define Materials_h 1

class G4Material;
class G4NistManager;

class Materials
{
  public:

    ~Materials();

    static Materials* GetInstance();

    G4Material* GetMaterial(G4String);

  private:

    Materials();

    void CreateMaterials();

  private:

    static Materials* instance;

    G4NistManager*     nistMan;

};

#endif /*Materials_h*/
