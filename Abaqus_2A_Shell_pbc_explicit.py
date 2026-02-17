class PBCContext(object):
    """
    PBC helper, written in a style close to the original script:
      - builds coord_to_label using rounded coords
      - anchor node sets (ANCH_ORIGIN, ANCH_X, ANCH_Y)
      - creates RPs (RPHX/RPHY/RPKX/RPKY) in the same positions as old code
      - creates equations named like: PBC_X_U1_000001, etc.
      - applies macro BCs on RPs (H,K)
    External usage remains:
        ctx = PBCContext(model, a, inst, dec=5)
        info = ctx.setup(step_name=STEP_NAME, H=H, K=K, make_anchors=True)
    """
    def __init__(self, model, a, inst, dec=5, tol=None):
        self.model = model
        self.a = a
        self.inst = inst
        self.dec = dec
        self.tol = (10.0 ** (-dec)) if tol is None else tol
        # cache for 1-node sets
        self.set_cache = {}
        # bounds + sizes
        self._detect_bounds()
        # rounded coordinate lookup
        self.coord_to_label = {}
        self._build_coord_to_label()

    # ----------------------------
    # Helpers (close to your style)
    # ----------------------------
    def is_close(self, a, b):
        return abs(a - b) <= self.tol

    def nrm(self, val):
        # round like your nrm(val, L) but now relative to bounds
        r = round(val, self.dec)
        if abs(r) < self.tol:
            r = 0.0
        return r

    # ----------------------------
    # bounds (xmin..xmax etc.)
    # ----------------------------
    def _detect_bounds(self):
        xs = [nd.coordinates[0] for nd in self.inst.nodes]
        ys = [nd.coordinates[1] for nd in self.inst.nodes]
        zs = [nd.coordinates[2] for nd in self.inst.nodes]
        self.xmin, self.xmax = min(xs), max(xs)
        self.ymin, self.ymax = min(ys), max(ys)
        self.zmin, self.zmax = min(zs), max(zs)
        self.Lx = self.xmax - self.xmin
        self.Ly = self.ymax - self.ymin
        self.Lz = self.zmax - self.zmin
        # centers (like X0_center, etc.)
        self.X0_center = self.xmin + 0.5 * self.Lx
        self.Y0_center = self.ymin + 0.5 * self.Ly
        self.Z0_center = self.zmin + 0.5 * self.Lz

    def _key_from_xyz(self, x, y, z):
        # Make a "0..L" style key by shifting to xmin/ymin/zmin then rounding
        xr = self.nrm(x - self.xmin)
        yr = self.nrm(y - self.ymin)
        zr = self.nrm(z - self.zmin)
        # snap to 0 or L like your old nrm(val, L)
        if abs(xr - 0.0) < self.tol: xr = 0.0
        if abs(xr - self.Lx) < self.tol: xr = self.Lx
        if abs(yr - 0.0) < self.tol: yr = 0.0
        if abs(yr - self.Ly) < self.tol: yr = self.Ly
        if abs(zr - 0.0) < self.tol: zr = 0.0
        if abs(zr - self.Lz) < self.tol: zr = self.Lz
        return (xr, yr, zr)

    def _build_coord_to_label(self):
        self.coord_to_label = {}
        for nd in self.inst.nodes:
            x, y, z = nd.coordinates
            key = self._key_from_xyz(x, y, z)
            self.coord_to_label[key] = nd.label

    # ----------------------------
    # node_set_name (same style)
    # ----------------------------
    def node_set_name(self, label):
        if label in self.set_cache:
            return self.set_cache[label]
        sname = 'N_%08d' % label
        if sname not in self.a.sets.keys():
            self.a.Set(name=sname, nodes=self.inst.nodes.sequenceFromLabels((label,)))
        self.set_cache[label] = sname
        return sname

    # ----------------------------
    # Anchors (same names as old)
    # ----------------------------
    def apply_minimal_anchors(self, bc_step='Initial'):
        # Use "0..L" keys
        def make_1node_set(name, key_0L):
            if key_0L not in self.coord_to_label:
                raise RuntimeError("Could not find node at key=%s" % (str(key_0L),))
            lab = self.coord_to_label[key_0L]
            if name not in self.a.sets.keys():
                self.a.Set(name=name, nodes=self.inst.nodes.sequenceFromLabels((lab,)))
            return name
        make_1node_set('ANCH_ORIGIN', (0.0, 0.0, 0.0))
        if 'BC_Origin' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(name='BC_Origin', createStepName=bc_step,
                                      region=self.a.sets['ANCH_ORIGIN'],
                                      u1=0.0, u2=0.0, u3=0.0)
        if (self.Lx, 0.0, 0.0) in self.coord_to_label:
            make_1node_set('ANCH_X', (self.Lx, 0.0, 0.0))
            if 'BC_AnchorX' not in self.model.boundaryConditions.keys():
                self.model.DisplacementBC(name='BC_AnchorX', createStepName=bc_step,
                                          region=self.a.sets['ANCH_X'],
                                          u2=0.0, u3=0.0)
        if (0.0, self.Ly, 0.0) in self.coord_to_label:
            make_1node_set('ANCH_Y', (0.0, self.Ly, 0.0))
            if 'BC_AnchorY' not in self.model.boundaryConditions.keys():
                self.model.DisplacementBC(name='BC_AnchorY', createStepName=bc_step,
                                          region=self.a.sets['ANCH_Y'],
                                          u3=0.0)

    # ----------------------------
    # RPs (same positions as old)
    # ----------------------------
    def ensure_macro_rps_legacy(self):
        """
        Create RPs at the same positions as the original script,
        but shifted by xmin/ymin/zmin so it works even if the RVE is translated.
        """
        def ensure_set(set_name, point_xyz):
            if set_name in self.a.sets.keys():
                return
            rp = self.a.ReferencePoint(point=point_xyz)
            self.a.Set(name=set_name, referencePoints=(self.a.referencePoints[rp.id],))
        # old: (1.5*Lx,0,0) etc in 0..L coords
        ensure_set('RPHX', (self.xmin + 1.5*self.Lx, self.ymin + 0.0,        self.zmin + 0.0))
        ensure_set('RPHY', (self.xmin + 0.0,        self.ymin + 1.5*self.Ly, self.zmin + 0.0))
        ensure_set('RPKX', (self.xmin + 1.5*self.Lx, self.ymin + 1.5*self.Ly, self.zmin + 0.0))
        ensure_set('RPKY', (self.xmin + 1.5*self.Lx, self.ymin + 1.5*self.Ly, self.zmin + 0.2*self.Lz))

    # ----------------------------
    # PBC equations (same naming)
    # ----------------------------
    def apply_pbc_equations_legacy(self,
                                   PBC_X_U1=True, PBC_X_U2=True, PBC_X_U3=True,
                                   PBC_Y_U1=True, PBC_Y_U2=True, PBC_Y_U3=True,
                                   exclude_edge_on_y=True):
        # Make sure RPs exist
        self.ensure_macro_rps_legacy()
        # --- X-plane: nodes at x = Lx (in 0..L coords)
        eqn_count_x_u1 = 0
        eqn_count_x_u2 = 0
        eqn_count_x_u3 = 0
        # bounding box for x ~ xmax in real coords
        x_nodes = self.inst.nodes.getByBoundingBox(self.xmax - self.tol,
                                                   self.ymin - self.tol,
                                                   self.zmin - self.tol,
                                                   self.xmax + self.tol,
                                                   self.ymax + self.tol,
                                                   self.zmax + self.tol)
        for nd in x_nodes:
            x, y, z = nd.coordinates
            # convert to 0..L coords for pairing
            xr, yr, zr = self._key_from_xyz(x, y, z)
            yrel = (self.ymin + yr) - self.Y0_center
            zrel = (self.zmin + zr) - self.Z0_center
            partner_key = (0.0, yr, zr)
            if partner_key not in self.coord_to_label:
                raise RuntimeError("Missing x-partner for key=%s" % (str(partner_key),))
            labA = nd.label
            labB = self.coord_to_label[partner_key]
            setA = self.node_set_name(labA)
            setB = self.node_set_name(labB)
            # X-U1
            if PBC_X_U1:
                terms = (
                    ( 1.0, setA, 1),
                    (-1.0, setB, 1),
                    (-self.Lx*1.0, 'RPHX', 1),
                    ( self.Lx*zrel, 'RPKX', 1),
                )
                eqn_count_x_u1 += 1
                self.model.Equation(name='PBC_X_U1_%06d' % eqn_count_x_u1, terms=terms)
            # X-U2
            if PBC_X_U2:
                terms = (
                    ( 1.0, setA, 2),
                    (-1.0, setB, 2),
                    (-self.Lx*1.0, 'RPHX', 2),
                    ( self.Lx*zrel, 'RPKX', 2),
                )
                eqn_count_x_u2 += 1
                self.model.Equation(name='PBC_X_U2_%06d' % eqn_count_x_u2, terms=terms)
            # X-U3
            if PBC_X_U3:
                terms = (
                    ( 1.0, setA, 3),
                    (-1.0, setB, 3),
                    (-self.Lx*1.0, 'RPHX', 3),
                    (-self.Lx*yrel, 'RPKX', 2),
                )
                eqn_count_x_u3 += 1
                self.model.Equation(name='PBC_X_U3_%06d' % eqn_count_x_u3, terms=terms)
        # --- Y-plane: nodes at y = Ly (in 0..L coords)
        eqn_count_y_u1 = 0
        eqn_count_y_u2 = 0
        eqn_count_y_u3 = 0
        y_nodes_all = self.inst.nodes.getByBoundingBox(self.xmin - self.tol,
                                                       self.ymax - self.tol,
                                                       self.zmin - self.tol,
                                                       self.xmax + self.tol,
                                                       self.ymax + self.tol,
                                                       self.zmax + self.tol)
        if exclude_edge_on_y:
            y_nodes = []
            for nd in y_nodes_all:
                xr, yr, zr = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
                if not self.is_close(xr, self.Lx):
                    y_nodes.append(nd)
        else:
            y_nodes = y_nodes_all
        for nd in y_nodes:
            x, y, z = nd.coordinates
            xr, yr, zr = self._key_from_xyz(x, y, z)
            xrel = (self.xmin + xr) - self.X0_center
            zrel = (self.zmin + zr) - self.Z0_center
            partner_key = (xr, 0.0, zr)
            if partner_key not in self.coord_to_label:
                raise RuntimeError("Missing y-partner for key=%s" % (str(partner_key),))
            labA = nd.label
            labB = self.coord_to_label[partner_key]
            setA = self.node_set_name(labA)
            setB = self.node_set_name(labB)
            # Y-U1
            if PBC_Y_U1:
                terms = (
                    ( 1.0, setA, 1),
                    (-1.0, setB, 1),
                    (-self.Ly*1.0, 'RPHY', 1),
                    (+self.Ly*zrel, 'RPKX', 2),
                )
                eqn_count_y_u1 += 1
                self.model.Equation(name='PBC_Y_U1_%06d' % eqn_count_y_u1, terms=terms)
            # Y-U2
            if PBC_Y_U2:
                terms = (
                    ( 1.0, setA, 2),
                    (-1.0, setB, 2),
                    (-self.Ly*1.0, 'RPHY', 2),
                    (+self.Ly*zrel, 'RPKY', 1),
                )
                eqn_count_y_u2 += 1
                self.model.Equation(name='PBC_Y_U2_%06d' % eqn_count_y_u2, terms=terms)
            # Y-U3
            if PBC_Y_U3:
                terms = (
                    ( 1.0, setA, 3),
                    (-1.0, setB, 3),
                    (-self.Ly*1.0, 'RPHY', 3),
                    (-self.Ly*xrel, 'RPKX', 2),
                )
                eqn_count_y_u3 += 1
                self.model.Equation(name='PBC_Y_U3_%06d' % eqn_count_y_u3, terms=terms)
        return {
            'RP_sets': {'RPHX': 'RPHX', 'RPHY': 'RPHY', 'RPKX': 'RPKX', 'RPKY': 'RPKY'},
            'counts': {
                'x_u1': eqn_count_x_u1, 'x_u2': eqn_count_x_u2, 'x_u3': eqn_count_x_u3,
                'y_u1': eqn_count_y_u1, 'y_u2': eqn_count_y_u2, 'y_u3': eqn_count_y_u3,
                'total': (eqn_count_x_u1 + eqn_count_x_u2 + eqn_count_x_u3 +
                          eqn_count_y_u1 + eqn_count_y_u2 + eqn_count_y_u3)
            }
        }

    # ----------------------------
    # Macro BCs (keep your logic)
    # ----------------------------
    def apply_macro_BCs(self, step_name, H, K, amplitude=None):
        # same as your original logic, but uses sets created here
        if 'BC_RPHX' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(
                name='BC_RPHX', createStepName=step_name,
                region=self.a.sets['RPHX'],
                u1=(H['H11'] - 1.0),
                u2=(H['H21'] - 0.0),
                u3=(H['H31'] - 0.0),
                amplitude=amplitude
            )
        if 'BC_RPHY' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(
                name='BC_RPHY', createStepName=step_name,
                region=self.a.sets['RPHY'],
                u1=(H['H12'] - 0.0),
                u2=(H['H22'] - 1.0),
                u3=(H['H32'] - 0.0),
                amplitude=amplitude
            )
        if 'BC_RPKX' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(
                name='BC_RPKX', createStepName=step_name,
                region=self.a.sets['RPKX'],
                u1=K['K11'], u2=K['K12'], u3=0.0,
                amplitude=amplitude
            )
        if 'BC_RPKY' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(
                name='BC_RPKY', createStepName=step_name,
                region=self.a.sets['RPKY'],
                u1=K['K22'], u2=0.0, u3=0.0,
                amplitude=amplitude
            )

    # ----------------------------
    # setup (same call signature)
    # ----------------------------
    def setup(self, step_name, H, K,
              make_anchors=True,
              anchors_step='Initial',
              exclude_edge_on_y=True,
              name_prefix='PBC',
              # switches exposed but optional (defaults match your old)
              PBC_X_U1=True, PBC_X_U2=True, PBC_X_U3=True,
              PBC_Y_U1=True, PBC_Y_U2=True, PBC_Y_U3=True):
        """
        One call: anchors + PBC equations + macro BCs.
        Keeps your external usage the same (H,K dicts).
        """
        if make_anchors:
            self.apply_minimal_anchors(bc_step=anchors_step)
        info = self.apply_pbc_equations_legacy(
            PBC_X_U1=PBC_X_U1, PBC_X_U2=PBC_X_U2, PBC_X_U3=PBC_X_U3,
            PBC_Y_U1=PBC_Y_U1, PBC_Y_U2=PBC_Y_U2, PBC_Y_U3=PBC_Y_U3,
            exclude_edge_on_y=exclude_edge_on_y
        )
        # Apply macro loads
        self.apply_macro_BCs(step_name=step_name, H=H, K=K, amplitude=None)
        return info



    # ==========================================================
    # Non-uniform nodal integration weights (dx*dy*dz)
    # Structured grid with nonuniform spacing in z (like your code)
    # ==========================================================
    def _build_axis_grids_0L(self):
        """
        Build sorted unique coordinate lists in 0..L coordinates (translation-safe),
        plus index maps for fast access.
        """
        # 0..L coords
        xvals = sorted(set(float(self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])[0])
                           for nd in self.inst.nodes))
        yvals = sorted(set(float(self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])[1])
                           for nd in self.inst.nodes))
        zvals = sorted(set(float(self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])[2])
                           for nd in self.inst.nodes))
        x_index = {x: i for i, x in enumerate(xvals)}
        y_index = {y: i for i, y in enumerate(yvals)}
        z_index = {z: i for i, z in enumerate(zvals)}
        return xvals, yvals, zvals, x_index, y_index, z_index

    def _local_delta(self, vals, idx):
        n = len(vals)
        if n < 2:
            return 1.0
        if idx == 0:
            return 0.5 * (vals[1] - vals[0])
        if idx == n - 1:
            return 0.5 * (vals[-1] - vals[-2])
        return 0.5 * (vals[idx + 1] - vals[idx - 1])

    def nodal_volume_weight_0L(self, x0L, y0L, z0L, grids=None):
        """
        dx*dy*dz weight at a node given (x,y,z) in 0..L coordinates.
        """
        if grids is None:
            grids = self._build_axis_grids_0L()
        xvals, yvals, zvals, x_index, y_index, z_index = grids
        ix = x_index[x0L]
        iy = y_index[y0L]
        iz = z_index[z0L]
        dx = self._local_delta(xvals, ix)
        dy = self._local_delta(yvals, iy)
        dz = self._local_delta(zvals, iz)
        return dx * dy * dz

    # ==========================================================
    # Pick an interior dependent node for moment constraint
    # ==========================================================
    def pick_interior_dep_node_for_moment(self, ztol=0.0):
        """
        Pick a TRUE interior node (not on boundary faces) and with nonzero
        coefficient wvol*(z - Z0_center).
        Uses bounding-box in real coords but evaluates coefficient in 0..L coords.
        """
        interior_nodes = self.inst.nodes.getByBoundingBox(
            self.xmin + self.tol, self.ymin + self.tol, self.zmin + self.tol,
            self.xmax - self.tol, self.ymax - self.tol, self.zmax - self.tol
        )
        if len(interior_nodes) == 0:
            raise RuntimeError("No interior nodes found. Refine mesh or reduce tol.")
        grids = self._build_axis_grids_0L()
        best_nd = None
        best_abscoef = -1.0
        for nd in interior_nodes:
            x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
            zrel = (self.zmin + z0L) - self.Z0_center  # real z - center
            if abs(zrel) <= max(ztol, 1e-30):
                continue
            wvol = self.nodal_volume_weight_0L(x0L, y0L, z0L, grids=grids)
            coef = wvol * zrel
            if abs(coef) > best_abscoef:
                best_abscoef = abs(coef)
                best_nd = nd
        if best_nd is None:
            raise RuntimeError("Interior nodes exist but all have ~zero moment coefficient. "
                               "Refine z or check Z0_center.")
        return best_nd

    # ==========================================================
    # GLOBAL zero first-moment constraint, weighted by dx*dy*dz
    # ==========================================================
    def add_zero_first_moment_weighted(self, eq_name, dof, ztol=0.0, echo=True):
        """
        Creates one Equation:
          sum_i wvol_i * (z_i - Z0_center) * U_dof(i) = 0
        Implemented by selecting ONE interior dependent node and moving its term first.
        """
        dep_node = self.pick_interior_dep_node_for_moment(ztol=ztol)
        dep_label = dep_node.label
        grids = self._build_axis_grids_0L()
        dep_coef = None
        terms_rest = []
        for nd in self.inst.nodes:
            x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
            zrel = (self.zmin + z0L) - self.Z0_center
            if abs(zrel) <= ztol:
                continue
            wvol = self.nodal_volume_weight_0L(x0L, y0L, z0L, grids=grids)
            coef = wvol * zrel
            if abs(coef) < 1e-30:
                continue
            if nd.label == dep_label:
                dep_coef = coef
            else:
                terms_rest.append((coef, self.node_set_name(nd.label), dof))
        if dep_coef is None or abs(dep_coef) < 1e-30:
            raise RuntimeError("Chosen dependent node has zero coefficient in %s." % eq_name)
        # dependent term MUST be first
        terms = [(dep_coef, self.node_set_name(dep_label), dof)] + terms_rest
        self.model.Equation(name=eq_name, terms=tuple(terms))
        if echo:
            print("Added:", eq_name,
                  " dependent node =", dep_label,
                  " dep_coef =", dep_coef,
                  " nTerms =", len(terms))




    # ==========================================================
    # Symmetric YZ-plane first-moment constraints (surface integral)
    # Implements the "image rule":
    #   - include also x=0 and x=Lx
    #   - if center plane exists, apply once
    #   - pairing tolerance = max(self.tol, 0.01 * min(dx_element))
    #
    # Continuous constraint:
    #   ∫_{x=xi} (z-z0) u_dof dA + ∫_{x=Lx-xi} (z-z0) u_dof dA = 0
    #
    # Discrete nodal quadrature:
    #   sum_nodes coef(node) * U_dof(node) = 0
    # with coef = (z-z0)*dA, dA = dy*dz from structured grid deltas.
    # ==========================================================

    def _x0L_of_node(self, nd):
        x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
        return float(x0L)

    def _y0Lz0L_of_node(self, nd):
        x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
        return float(y0L), float(z0L)

    def xplanes_0L(self):
        """
        Sorted unique x-planes in 0..L coordinates (already rounded/snapped by _key_from_xyz).
        """
        return sorted(set(self._x0L_of_node(nd) for nd in self.inst.nodes))

    def _dx_min_xplanes(self, xvals=None):
        """
        Minimum spacing between x-planes (acts as dx_element).
        """
        if xvals is None:
            xvals = self.xplanes_0L()
        if len(xvals) < 2:
            return 1.0
        dxs = [xvals[i+1] - xvals[i] for i in range(len(xvals) - 1)]
        return min(dxs) if dxs else 1.0

    def _plane_pair_tol(self, xvals=None):
        """
        Pairing tolerance from image rule:
            tol_pair = 0.01 * min(dx_element)
        but never smaller than self.tol (your rounding tolerance).
        """
        dxmin = self._dx_min_xplanes(xvals=xvals)
        return max(self.tol, 0.01 * dxmin)

    def _find_nearest_plane(self, x_target, xvals, tolx):
        """
        Return the plane value in xvals that is closest to x_target, if within tolx.
        """
        xt = float(x_target)
        best = None
        best_err = 1e100
        for x in xvals:
            err = abs(float(x) - xt)
            if err < best_err:
                best_err = err
                best = float(x)
        if best is None or best_err > tolx:
            return None
        return best

    def collect_xplane_nodes_0L(self, xface0L):
        """
        Collect nodes on x = xface0L using tolerance-based match (robust).
        """
        xf = float(xface0L)
        xvals = self.xplanes_0L()
        tolx = self._plane_pair_tol(xvals=xvals)

        nodes = []
        for nd in self.inst.nodes:
            xr = self._x0L_of_node(nd)
            if abs(xr - xf) <= tolx:
                nodes.append(nd)
        return nodes

    def symmetric_xplane_0L(self, xface0L):
        """
        For a plane xi, find its symmetric partner xj_target = Lx - xi,
        then match to the nearest existing plane within tol_pair.
        """
        xi = float(xface0L)
        xvals = self.xplanes_0L()
        tolx = self._plane_pair_tol(xvals=xvals)

        xj_target = float(self.Lx - xi)  # do NOT round here
        xj = self._find_nearest_plane(xj_target, xvals=xvals, tolx=tolx)
        if xj is None:
            raise RuntimeError(
                "Symmetric plane not found: xi=%g -> xj_target=%g, tol_pair=%g, xplanes=%s"
                % (xi, xj_target, tolx, xvals)
            )
        return xj

    def nodal_area_weight_yz_0L(self, y0L, z0L, grids=None):
        """
        dy*dz nodal area weight on a YZ plane, using structured-grid deltas.
        """
        if grids is None:
            grids = self._build_axis_grids_0L()
        xvals, yvals, zvals, x_index, y_index, z_index = grids

        iy = y_index[float(y0L)]
        iz = z_index[float(z0L)]
        dy = self._local_delta(yvals, iy)
        dz = self._local_delta(zvals, iz)
        return dy * dz

    def pick_interior_dep_node_for_yz_planes(self, nodes, ztol=0.0, require_nonzero_coef=True):
        """
        Dependent-node picker for YZ-plane moment equation.
    
        New rule:
          Pick the interior (in Y and Z) node with minimum y, then minimum z
          (in 0..L coordinates), i.e. lexicographic min (y0L, z0L).
    
        Interior definition:
          y not on {0, Ly} and z not on {0, Lz} (using self.is_close).
    
        Optionally still require nonzero coefficient dA*(z-z0) so the dependent
        node actually appears in the moment equation.
        """
        grids = self._build_axis_grids_0L()
    
        def is_interior_yz(y0L, z0L):
            if self.is_close(y0L, 0.0) or self.is_close(y0L, self.Ly):
                return False
            if self.is_close(z0L, 0.0) or self.is_close(z0L, self.Lz):
                return False
            return True
    
        best_nd = None
        best_key = None  # (y, z)
    
        for nd in nodes:
            y0L, z0L = self._y0Lz0L_of_node(nd)
    
            # Only interior nodes in Y and Z
            if not is_interior_yz(y0L, z0L):
                continue
    
            if require_nonzero_coef:
                zrel = (self.zmin + z0L) - self.Z0_center
                if abs(zrel) <= max(ztol, 1e-30):
                    continue
                dA = self.nodal_area_weight_yz_0L(y0L, z0L, grids=grids)
                coef = dA * zrel
                if abs(coef) < 1e-30:
                    continue
    
            key = (float(y0L), float(z0L))
            if best_key is None or key < best_key:
                best_key = key
                best_nd = nd
    
        if best_nd is None:
            raise RuntimeError(
                "No interior dependent node found with min(y),min(z). "
                "Check that the plane has interior nodes (0<y<Ly, 0<z<Lz) "
                "and that ztol is not filtering everything."
            )
    
        return best_nd

    
        # PASS 1: your preferred rule
        best_nd = pick(require_interior_yz=True, require_x_thresh=True)
        if best_nd is not None:
            return best_nd
    
        # PASS 2: relax x restriction (still interior in Y,Z)
        best_nd = pick(require_interior_yz=True, require_x_thresh=False)
        if best_nd is not None:
            return best_nd
    
        # PASS 3: last resort (drop interior YZ requirement)
        best_nd = pick(require_interior_yz=False, require_x_thresh=False)
        if best_nd is not None:
            return best_nd
    
        raise RuntimeError(
            "No dependent node found: all coefficients filtered out (ztol too large?) "
            "or plane has no usable nodes."
        )



    def add_symmetric_first_moment_yz_planes(self, eq_name, xi_0L, dof=1, ztol=0.0, echo=True):
        """
        Build ONE equation for the pair (xi, xj=Lx-xi).
        If xi is the center plane (xi==xj), apply only once.
        Includes boundaries if xi is boundary (this method allows it).
        """
        xvals = self.xplanes_0L()
        tolx = self._plane_pair_tol(xvals=xvals)

        # snap xi to nearest existing plane
        xi = self._find_nearest_plane(float(xi_0L), xvals=xvals, tolx=tolx)
        if xi is None:
            raise RuntimeError("xi=%g is not an existing x-plane within tol_pair=%g." % (float(xi_0L), tolx))

        xj = self.symmetric_xplane_0L(xi)
        same_plane = abs(xi - xj) <= tolx

        nodes_i = self.collect_xplane_nodes_0L(xi)
        nodes_j = [] if same_plane else self.collect_xplane_nodes_0L(xj)

        grids = self._build_axis_grids_0L()
        coef_by_label = {}

        def accumulate(nodes):
            for nd in nodes:
                y0L, z0L = self._y0Lz0L_of_node(nd)
                zrel = (self.zmin + z0L) - self.Z0_center
                if abs(zrel) <= ztol:
                    continue
                dA = self.nodal_area_weight_yz_0L(y0L, z0L, grids=grids)
                coef = dA * zrel
                if abs(coef) < 1e-30:
                    continue
                lab = nd.label
                coef_by_label[lab] = coef_by_label.get(lab, 0.0) + coef

        accumulate(nodes_i)
        accumulate(nodes_j)

        if not coef_by_label:
            raise RuntimeError("All coefficients ~0 on xi=%g,xj=%g. Reduce ztol or refine mesh." % (xi, xj))

        # dependent node must come from union of nodes on both planes
        union_nodes = nodes_i + nodes_j
        dep_node = self.pick_interior_dep_node_for_yz_planes(union_nodes, ztol=ztol)
        dep_label = dep_node.label

        if dep_label not in coef_by_label:
            raise RuntimeError("Dependent node got filtered out (ztol too large?).")

        dep_coef = coef_by_label[dep_label]
        if abs(dep_coef) < 1e-30:
            raise RuntimeError("Dependent node coefficient is ~0 in %s." % eq_name)

        # dependent term first
        terms = [(dep_coef, self.node_set_name(dep_label), dof)]
        for lab, coef in coef_by_label.items():
            if lab == dep_label:
                continue
            terms.append((coef, self.node_set_name(lab), dof))

        self.model.Equation(name=eq_name, terms=tuple(terms))

        if echo:
            print("Added:", eq_name,
                  " xi=", xi, " xj=", xj, " same_plane=", same_plane,
                  " dep_label=", dep_label, " nTerms=", len(terms))

    def add_symmetric_first_moment_yz_all_pairs(self,
                                                eq_prefix="SYM_M1_UX_YZ",
                                                dof=1,
                                                ztol=0.0,
                                                echo=True):
        """
        Add ONE equation for EACH UNIQUE symmetric pair (xi, xj=Lx-xi),
        INCLUDING boundary planes x=0 and x=Lx.

        Expected number of equations:
            ceil(nPlanes/2) = (nPlanes + 1)//2
        """
        xvals = self.xplanes_0L()
        tolx = self._plane_pair_tol(xvals=xvals)

        used = set()
        k = 0
        skipped = 0

        for xi in xvals:
                # EXCLUDE boundary planes x=0 and x=Lx
            if self.is_close(float(xi), 0.0) or self.is_close(float(xi), self.Lx):
                continue
            try:
                xj = self.symmetric_xplane_0L(xi)
            except Exception:
                skipped += 1
                if echo:
                    print("Skip xi=%g (no symmetric partner within tol_pair=%g)" % (float(xi), tolx))
                continue

            a = min(float(xi), float(xj))
            b = max(float(xi), float(xj))
            # (optional but recommended) also skip if partner is boundary
            if self.is_close(a, 0.0) or self.is_close(b, self.Lx):
                continue
            
            # canonical bin key to avoid duplicates
            key = (int(round(a / tolx)), int(round(b / tolx)))
            if key in used:
                continue
            used.add(key)

            k += 1
            eq_name = "%s_%03d" % (eq_prefix, k)
            self.add_symmetric_first_moment_yz_planes(eq_name=eq_name,
                                                     xi_0L=a,
                                                     dof=dof,
                                                     ztol=ztol,
                                                     echo=echo)

        if echo:
            nPlanes = len(xvals)
            expected = (nPlanes + 1) // 2
            print("Added", k, "equation(s). Expected", expected,
                  "| nPlanes =", nPlanes,
                  "| skipped planes =", skipped,
                  "| tol_pair =", tolx)

        return k



    def _y0L_of_node(self, nd):
        x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
        return float(y0L)
    
    def _x0Lz0L_of_node(self, nd):
        x0L, y0L, z0L = self._key_from_xyz(nd.coordinates[0], nd.coordinates[1], nd.coordinates[2])
        return float(x0L), float(z0L)
    
    def yplanes_0L(self):
        """Sorted unique y-planes in 0..L coordinates."""
        return sorted(set(self._y0L_of_node(nd) for nd in self.inst.nodes))
    
    def _dy_min_yplanes(self, yvals=None):
        """Minimum spacing between y-planes (acts as dy_element)."""
        if yvals is None:
            yvals = self.yplanes_0L()
        if len(yvals) < 2:
            return 1.0
        dys = [yvals[i+1] - yvals[i] for i in range(len(yvals) - 1)]
        return min(dys) if dys else 1.0
    
    def _plane_pair_tol_y(self, yvals=None):
        """Pairing tolerance for y-plane symmetry (same rule as x): max(self.tol, 0.01*min(dy))."""
        dymin = self._dy_min_yplanes(yvals=yvals)
        return max(self.tol, 0.01 * dymin)
    
    def collect_yplane_nodes_0L(self, yface0L):
        """Collect nodes on y = yface0L using tolerance-based match."""
        yf = float(yface0L)
        yvals = self.yplanes_0L()
        toly = self._plane_pair_tol_y(yvals=yvals)
    
        nodes = []
        for nd in self.inst.nodes:
            yr = self._y0L_of_node(nd)
            if abs(yr - yf) <= toly:
                nodes.append(nd)
        return nodes
    
    def symmetric_yplane_0L(self, yface0L):
        """For a plane yi, find symmetric partner yj = Ly - yi, matched to nearest existing plane."""
        yi = float(yface0L)
        yvals = self.yplanes_0L()
        toly = self._plane_pair_tol_y(yvals=yvals)
    
        yj_target = float(self.Ly - yi)
        yj = self._find_nearest_plane(yj_target, xvals=yvals, tolx=toly)  # reuse your _find_nearest_plane
        if yj is None:
            raise RuntimeError(
                "Symmetric y-plane not found: yi=%g -> yj_target=%g, tol_pair=%g, yplanes=%s"
                % (yi, yj_target, toly, yvals)
            )
        return yj
    
    
    
    def nodal_area_weight_xz_0L(self, x0L, z0L, grids=None):
        """
        dx*dz nodal area weight on an XZ plane (y = const),
        using structured-grid deltas.
        """
        if grids is None:
            grids = self._build_axis_grids_0L()
        xvals, yvals, zvals, x_index, y_index, z_index = grids
    
        ix = x_index[float(x0L)]
        iz = z_index[float(z0L)]
        dx = self._local_delta(xvals, ix)
        dz = self._local_delta(zvals, iz)
        return dx * dz
    
    
    
    def pick_interior_dep_node_for_xz_planes(self, nodes, ztol=0.0, require_nonzero_coef=True):
        """
        Dependent-node picker for XZ-plane moment equation (y=const planes).
    
        New rule:
          Pick the interior (in X and Z) node with maximum x, then maximum z
          (in 0..L coordinates), i.e. lexicographic max (x0L, z0L).
    
        Interior definition:
          x not on {0, Lx} and z not on {0, Lz} (using self.is_close).
    
        Optionally require nonzero coefficient dA*(z-z0) (dA=dx*dz on XZ plane).
        This matches "u2*z integration" style constraints.
        """
        grids = self._build_axis_grids_0L()
    
        def is_interior_xz(x0L, z0L):
            if self.is_close(x0L, 0.0) or self.is_close(x0L, self.Lx):
                return False
            if self.is_close(z0L, 0.0) or self.is_close(z0L, self.Lz):
                return False
            return True
    
        best_nd = None
        best_key = None  # (x, z)
    
        for nd in nodes:
            x0L, z0L = self._x0Lz0L_of_node(nd)
    
            if not is_interior_xz(x0L, z0L):
                continue
    
            if require_nonzero_coef:
                zrel = (self.zmin + z0L) - self.Z0_center
                if abs(zrel) <= max(ztol, 1e-30):
                    continue
                dA = self.nodal_area_weight_xz_0L(x0L, z0L, grids=grids)
                coef = dA * zrel
                if abs(coef) < 1e-30:
                    continue
    
            key = (float(x0L), float(z0L))
            if best_key is None or key > best_key:
                best_key = key
                best_nd = nd
    
        if best_nd is None:
            raise RuntimeError(
                "No interior dependent node found with max(x),max(z). "
                "Check that the plane has interior nodes (0<x<Lx, 0<z<Lz) "
                "and that ztol is not filtering everything."
            )
    
        return best_nd




    def add_symmetric_first_moment_xz_planes(self, eq_name, yi_0L, dof=2, ztol=0.0, echo=True):
        """
        Mirror of add_symmetric_first_moment_yz_planes, but for XZ planes (y=const).
    
        Continuous idea:
          ∫_{y=yi} (z-z0) u2 dA + ∫_{y=Ly-yi} (z-z0) u2 dA = 0
    
        Discrete:
          sum_nodes coef(node)*U_dof(node) = 0
          coef = (z-z0)*dA,  dA = dx*dz
        """
        yvals = self.yplanes_0L()
        toly = self._plane_pair_tol_y(yvals=yvals)
    
        # snap yi to nearest existing y-plane
        yi = self._find_nearest_plane(float(yi_0L), xvals=yvals, tolx=toly)
        if yi is None:
            raise RuntimeError("yi=%g is not an existing y-plane within tol_pair=%g." % (float(yi_0L), toly))
    
        yj = self.symmetric_yplane_0L(yi)
        same_plane = abs(float(yi) - float(yj)) <= toly
    
        nodes_i = self.collect_yplane_nodes_0L(yi)
        nodes_j = [] if same_plane else self.collect_yplane_nodes_0L(yj)
    
        grids = self._build_axis_grids_0L()
        coef_by_label = {}
    
        def accumulate(nodes):
            for nd in nodes:
                x0L, z0L = self._x0Lz0L_of_node(nd)
                zrel = (self.zmin + z0L) - self.Z0_center
                if abs(zrel) <= ztol:
                    continue
                dA = self.nodal_area_weight_xz_0L(x0L, z0L, grids=grids)
                coef = dA * zrel
                if abs(coef) < 1e-30:
                    continue
                lab = nd.label
                coef_by_label[lab] = coef_by_label.get(lab, 0.0) + coef
    
        accumulate(nodes_i)
        accumulate(nodes_j)
    
        if not coef_by_label:
            raise RuntimeError("All coefficients ~0 on yi=%g,yj=%g. Reduce ztol or refine mesh." % (yi, yj))
    
        # dependent node from union; deterministic max(x),max(z) interior
        union_nodes = nodes_i + nodes_j
        dep_node = self.pick_interior_dep_node_for_xz_planes(union_nodes, ztol=ztol)
        dep_label = dep_node.label
    
        if dep_label not in coef_by_label:
            raise RuntimeError("Dependent node got filtered out (ztol too large?).")
    
        dep_coef = coef_by_label[dep_label]
        if abs(dep_coef) < 1e-30:
            raise RuntimeError("Dependent node coefficient is ~0 in %s." % eq_name)
    
        # dependent term first
        terms = [(dep_coef, self.node_set_name(dep_label), dof)]
        for lab, coef in coef_by_label.items():
            if lab == dep_label:
                continue
            terms.append((coef, self.node_set_name(lab), dof))
    
        self.model.Equation(name=eq_name, terms=tuple(terms))
    
        if echo:
            print("Added:", eq_name,
                  " yi=", yi, " yj=", yj, " same_plane=", same_plane,
                  " dep_label=", dep_label, " nTerms=", len(terms))
    
    
    def add_symmetric_first_moment_xz_all_pairs(self,
                                                eq_prefix="SYM_M1_U2_XZ",
                                                dof=2,
                                                ztol=0.0,
                                                echo=True,
                                                exclude_boundary=True):
        """
        Mirror of add_symmetric_first_moment_yz_all_pairs, but for y-planes.
    
        If exclude_boundary=True: skip y=0 and y=Ly (matches your yz_all_pairs style).
        """
        yvals = self.yplanes_0L()
        toly = self._plane_pair_tol_y(yvals=yvals)
    
        used = set()
        k = 0
        skipped = 0
    
        for yi in yvals:
            if exclude_boundary and (self.is_close(float(yi), 0.0) or self.is_close(float(yi), self.Ly)):
                continue
    
            try:
                yj = self.symmetric_yplane_0L(yi)
            except Exception:
                skipped += 1
                if echo:
                    print("Skip yi=%g (no symmetric partner within tol_pair=%g)" % (float(yi), toly))
                continue
    
            a = min(float(yi), float(yj))
            b = max(float(yi), float(yj))
    
            if exclude_boundary and (self.is_close(a, 0.0) or self.is_close(b, self.Ly)):
                continue
    
            # canonical key to avoid duplicates
            key = (int(round(a / toly)), int(round(b / toly)))
            if key in used:
                continue
            used.add(key)
    
            k += 1
            eq_name = "%s_%03d" % (eq_prefix, k)
            self.add_symmetric_first_moment_xz_planes(eq_name=eq_name,
                                                     yi_0L=a,
                                                     dof=dof,
                                                     ztol=ztol,
                                                     echo=echo)
    
        if echo:
            nPlanes = len(yvals)
            expected = (nPlanes + 1) // 2
            print("Added", k, "equation(s). Expected", expected,
                  "| nPlanes =", nPlanes,
                  "| skipped planes =", skipped,
                  "| tol_pair =", toly)
    
        return k
    
