    
    std::shared_ptr<MeshFunction<size_t>> boundary_markers =
        std::make_shared<MeshFunction<size_t>>(mesh, 2);
    boundary_markers->set_all(9999);
    
    // setup boundary objects
    std::shared_ptr<BoundaryX0> bx0 = std::make_shared<BoundaryX0>();
    bx0->setLwh(l,w,h);
    std::shared_ptr<BoundaryX1> bx1 = std::make_shared<BoundaryX1>();
    bx1->setLwh(l,w,h);
    std::shared_ptr<BoundaryY0> by0 = std::make_shared<BoundaryY0>();
    by0->setLwh(l,w,h);
    std::shared_ptr<BoundaryY1> by1 = std::make_shared<BoundaryY1>();
    by1->setLwh(l,w,h);
    std::shared_ptr<BoundaryZ0> bz0 = std::make_shared<BoundaryZ0>();
    bz0->setLwh(l,w,h);
    std::shared_ptr<BoundaryZ1> bz1 = std::make_shared<BoundaryZ1>();
    bz1->setLwh(l,w,h);

    bx0->mark(*boundary_markers, 0);
    bx1->mark(*boundary_markers, 1);
    by0->mark(*boundary_markers, 2);
    by1->mark(*boundary_markers, 3);
    bz0->mark(*boundary_markers, 4);
    bz1->mark(*boundary_markers, 5);
