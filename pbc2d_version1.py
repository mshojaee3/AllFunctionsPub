# 2D_pbc.py
from abaqusConstants import *
from regionToolset import Region

class PBC2DContext(object):
    """
    2D Periodic Boundary Condition helper (RVE rectangle) in the same *style*
    as shell_pbc.py:
      - context object storing model/assembly/instance
      - rounded coord_to_label lookup
      - anchor node sets
      - one-node set cache (N_########)
      - setup() as one-call entry point
      - creates RPX/RPY + BCs on them (macro jump)
      - creates equations for:
          * Right-Left pairing (exclude corners)
          * Top-Bottom pairing (exclude corners)
          * Corner periodicity equations
    """

    def __init__(self, model, a, inst, dec=8, tol=None):
        self.model = model
        self.a = a
        self.inst = inst
        self.dec = dec
        self.tol = (10.0 ** (-dec)) if tol is None else tol

        self.set_cache = {}
        self.coord_to_label = {}

        self._detect_bounds()
        self._build_coord_to_label()

    # ----------------------------
    # geometry bounds
    # ----------------------------
    def _detect_bounds(self):
        xs = [nd.coordinates[0] for nd in self.inst.nodes]
        ys = [nd.coordinates[1] for nd in self.inst.nodes]
        self.xmin, self.xmax = min(xs), max(xs)
        self.ymin, self.ymax = min(ys), max(ys)
        self.Lx = self.xmax - self.xmin
        self.Ly = self.ymax - self.ymin

    # ----------------------------
    # rounding + keys
    # ----------------------------
    def nrm(self, v):
        r = round(v, self.dec)
        if abs(r) < self.tol:
            r = 0.0
        return r

    def key_xy_0L(self, x, y):
        # convert to 0..L coordinates, then round/snap
        xr = self.nrm(x - self.xmin)
        yr = self.nrm(y - self.ymin)

        # snap to {0, L}
        if abs(xr - 0.0) < self.tol: xr = 0.0
        if abs(xr - self.Lx) < self.tol: xr = self.Lx
        if abs(yr - 0.0) < self.tol: yr = 0.0
        if abs(yr - self.Ly) < self.tol: yr = self.Ly
        return (xr, yr)

    def _build_coord_to_label(self):
        self.coord_to_label = {}
        for nd in self.inst.nodes:
            x, y, z = nd.coordinates
            self.coord_to_label[self.key_xy_0L(x, y)] = nd.label

    # ----------------------------
    # node set cache (same style)
    # ----------------------------
    def node_set_name(self, label, prefix='N_'):
        if label in self.set_cache:
            return self.set_cache[label]
        sname = '%s%08d' % (prefix, label)
        if sname not in self.a.sets.keys():
            self.a.Set(name=sname, nodes=self.inst.nodes.sequenceFromLabels((label,)))
        self.set_cache[label] = sname
        return sname

    def make_1node_set_by_xy0L(self, set_name, x0L, y0L):
        key = (float(x0L), float(y0L))
        if key not in self.coord_to_label:
            raise RuntimeError("Node not found at (x0L,y0L)=%s. Try larger tol/dec." % (str(key),))
        lab = self.coord_to_label[key]
        if set_name not in self.a.sets.keys():
            self.a.Set(name=set_name, nodes=self.inst.nodes.sequenceFromLabels((lab,)))
        return set_name

    # ----------------------------
    # anchors (remove RBMs)
    # ----------------------------
    def apply_minimal_anchors(self, anchors_step='Initial',
                              name_origin='ANCH_00', name_x='ANCH_10'):
        """
        Fix origin: u1=u2=0
        Fix (Lx,0): u2=0  (prevents rigid rotation)
        """
        self.make_1node_set_by_xy0L(name_origin, 0.0, 0.0)
        if 'BC_00' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(name='BC_00', createStepName=anchors_step,
                                      region=self.a.sets[name_origin],
                                      u1=0.0, u2=0.0, ur3=UNSET)

        self.make_1node_set_by_xy0L(name_x, self.Lx, 0.0)
        if 'BC_10' not in self.model.boundaryConditions.keys():
            self.model.DisplacementBC(name='BC_10', createStepName=anchors_step,
                                      region=self.a.sets[name_x],
                                      u1=UNSET, u2=0.0, ur3=UNSET)

    # ----------------------------
    # RPs + macro BCs
    # ----------------------------
    def ensure_rps(self, rpx_name='RPX', rpy_name='RPY'):
        """
        Place RPs using the same idea as your script:
          RPX at (xmin + 1.5*Lx, ymin, 0)
          RPY at (xmin, ymin + 1.5*Ly, 0)
        """
        if rpx_name not in self.a.sets.keys():
            rpX = self.a.ReferencePoint(point=(self.xmin + 1.5*self.Lx, self.ymin + 0.0, 0.0))
            self.a.Set(name=rpx_name, referencePoints=(self.a.referencePoints[rpX.id],))
        if rpy_name not in self.a.sets.keys():
            rpY = self.a.ReferencePoint(point=(self.xmin + 0.0, self.ymin + 1.5*self.Ly, 0.0))
            self.a.Set(name=rpy_name, referencePoints=(self.a.referencePoints[rpY.id],))

    def apply_macro_BCs(self, step_name, H11, H12, H21, H22,
                        rpx_name='RPX', rpy_name='RPY',
                        bc_rpx='BC_RPX', bc_rpy='BC_RPY'):
        """
        Macro jump:
          u(RPX) = [H11*Lx, H21*Lx]
          u(RPY) = [H12*Ly, H22*Ly]
        """
        self.ensure_rps(rpx_name=rpx_name, rpy_name=rpy_name)

        if bc_rpx in self.model.boundaryConditions.keys():
            del self.model.boundaryConditions[bc_rpx]
        if bc_rpy in self.model.boundaryConditions.keys():
            del self.model.boundaryConditions[bc_rpy]

        self.model.DisplacementBC(name=bc_rpx, createStepName=step_name,
                                  region=self.a.sets[rpx_name],
                                  u1=H11*self.Lx, u2=H21*self.Lx, ur3=UNSET)
        self.model.DisplacementBC(name=bc_rpy, createStepName=step_name,
                                  region=self.a.sets[rpy_name],
                                  u1=H12*self.Ly, u2=H22*self.Ly, ur3=UNSET)

    # ----------------------------
    # build node lists on boundaries
    # ----------------------------
    def _collect_boundary_nodes(self, exclude_corners=True):
        right_nodes = []
        top_nodes = []

        for nd in self.inst.nodes:
            x, y, z = nd.coordinates
            x0L, y0L = self.key_xy_0L(x, y)

            on_right = abs(x0L - self.Lx) <= self.tol
            on_left  = abs(x0L - 0.0)    <= self.tol
            on_top   = abs(y0L - self.Ly) <= self.tol
            on_bot   = abs(y0L - 0.0)     <= self.tol

            if exclude_corners:
                # for LR pairing: skip y=0 and y=Ly
                if on_right and (y0L > self.tol) and (y0L < self.Ly - self.tol):
                    right_nodes.append((nd.label, x0L, y0L))
                # for TB pairing: skip x=0 and x=Lx
                if on_top and (x0L > self.tol) and (x0L < self.Lx - self.tol):
                    top_nodes.append((nd.label, x0L, y0L))
            else:
                if on_right and not on_left:
                    right_nodes.append((nd.label, x0L, y0L))
                if on_top and not on_bot:
                    top_nodes.append((nd.label, x0L, y0L))

        return right_nodes, top_nodes

    # ----------------------------
    # equations (2D style)
    # ----------------------------
    def apply_pbc_equations(self, rpx_name='RPX', rpy_name='RPY',
                            exclude_corners=True,
                            eq_prefix='PBC',
                            dofs=(1, 2)):
        """
        Create:
          uR - uL - uRPX = 0   for right/left pairs (same y)
          uT - uB - uRPY = 0   for top/bottom pairs (same x)
        Plus corner periodicity constraints (always).
        """
        self.ensure_rps(rpx_name=rpx_name, rpy_name=rpy_name)

        right_nodes, top_nodes = self._collect_boundary_nodes(exclude_corners=exclude_corners)

        # Right-Left pairs
        eqn = 0
        for labR, xR0L, yR0L in right_nodes:
            partner_key = (0.0, yR0L)
            if partner_key not in self.coord_to_label:
                raise RuntimeError("Missing LEFT partner for y0L=%g" % yR0L)
            labL = self.coord_to_label[partner_key]
            setR = self.node_set_name(labR)
            setL = self.node_set_name(labL)

            for dof in dofs:
                eqn += 1
                self.model.Equation(
                    name='%s_LR_%06d' % (eq_prefix, eqn),
                    terms=(( 1.0, setR, dof),
                           (-1.0, setL, dof),
                           (-1.0, rpx_name, dof))
                )

        # Top-Bottom pairs
        for labT, xT0L, yT0L in top_nodes:
            partner_key = (xT0L, 0.0)
            if partner_key not in self.coord_to_label:
                raise RuntimeError("Missing BOTTOM partner for x0L=%g" % xT0L)
            labB = self.coord_to_label[partner_key]
            setT = self.node_set_name(labT)
            setB = self.node_set_name(labB)

            for dof in dofs:
                eqn += 1
                self.model.Equation(
                    name='%s_TB_%06d' % (eq_prefix, eqn),
                    terms=(( 1.0, setT, dof),
                           (-1.0, setB, dof),
                           (-1.0, rpy_name, dof))
                )

        # Corners (one node each)
        self.make_1node_set_by_xy0L('P1', 0.0, 0.0)
        self.make_1node_set_by_xy0L('P2', self.Lx, 0.0)
        self.make_1node_set_by_xy0L('P4', 0.0, self.Ly)
        self.make_1node_set_by_xy0L('P3', self.Lx, self.Ly)

        # Corner periodicity:
        # P2 - P1 = RPX
        # P4 - P1 = RPY
        # P3 - P2 = RPY
        for dof in dofs:
            self.model.Equation(
                name='%s_P2P1_%d' % (eq_prefix, dof),
                terms=((+1.0, 'P2', dof),
                       (-1.0, 'P1', dof),
                       (-1.0, rpx_name, dof))
            )
            self.model.Equation(
                name='%s_P4P1_%d' % (eq_prefix, dof),
                terms=((+1.0, 'P4', dof),
                       (-1.0, 'P1', dof),
                       (-1.0, rpy_name, dof))
            )
            self.model.Equation(
                name='%s_P3P2_%d' % (eq_prefix, dof),
                terms=((+1.0, 'P3', dof),
                       (-1.0, 'P2', dof),
                       (-1.0, rpy_name, dof))
            )

        return {'eq_count_total': eqn, 'n_right': len(right_nodes), 'n_top': len(top_nodes)}

    # ----------------------------
    # setup (one-call API)
    # ----------------------------
    def setup(self,
              step_name,
              H11, H12, H21, H22,
              make_anchors=True,
              anchors_step='Initial',
              rpx_name='RPX',
              rpy_name='RPY',
              exclude_corners=True,
              eq_prefix='PBC'):
        """
        One call:
          - anchors
          - RPs + macro BCs
          - PBC equations + corner equations
        """
        if make_anchors:
            self.apply_minimal_anchors(anchors_step=anchors_step)

        self.apply_macro_BCs(step_name=step_name,
                             H11=H11, H12=H12, H21=H21, H22=H22,
                             rpx_name=rpx_name, rpy_name=rpy_name)

        info = self.apply_pbc_equations(rpx_name=rpx_name, rpy_name=rpy_name,
                                        exclude_corners=exclude_corners,
                                        eq_prefix=eq_prefix,
                                        dofs=(1, 2))

        print("2D PBC: created %d equations | right nodes=%d | top nodes=%d"
              % (info['eq_count_total'], info['n_right'], info['n_top']))
        return info
