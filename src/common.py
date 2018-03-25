from enum import Enum


class FilterType(Enum):
    """
    Filter type (filter densities, or gradients only)
    """
    NoFilter = 0
    Density = 1
    Sensitivity = 2


class InterpolationType(Enum):
    """
    Material interpolation scheme: classic SIMP, or Pedersen (for self-weight problems)
    """
    SIMP = 1
    Pedersen = 2


class ProblemType(Enum):
    """
    Problem type. Minimize appearance only, minimize compliance only, or minimize
    appearance with a compliance constraint.
    """
    Appearance = 1
    Compliance = 2
    AppearanceWithMaxCompliance = 3
    ###Added by Antoine Hoffmann EPFL 2018
    ComplianceWithSymmetry = 4
    AppearanceWithMaxComplianceAndSymmetry = 5
    
    def involves_symmetry(self):
		"""
		Returns true iff the given problem type has symmetry
		"""
		
		return self in (ProblemType.ComplianceWithSymmetry
		                ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)
	######
    def involves_appearance(self):
        """
        Returns true iff the given problem type requires the appearance evaluation.
        """
        return self in (ProblemType.Appearance, ProblemType.AppearanceWithMaxCompliance
                        ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)

    def involves_compliance(self):
        """
        Returns true iff the given problem type requires the appearance evaluation.
        """
        return self in (ProblemType.Compliance, ProblemType.AppearanceWithMaxCompliance
                       ,ProblemType.ComplianceWithSymmetry
                       ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)

    def involves_volume(self):
        """
        Returns true iff the given problem type has a volume constraint.
        """
        return self in (ProblemType.Compliance, ProblemType.AppearanceWithMaxCompliance
                       ,ProblemType.ComplianceWithSymmetry
                       ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)

    def has_compliance_constraint(self):
        """
        Returns true iff the given problem has a constraint on the compliance.
        """
        return self in (ProblemType.AppearanceWithMaxCompliance 
                       ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)

    def has_volume_constraint(self):
        """
        Returns true iff the given problem has a constraint on the volume.
        """
        return self in (ProblemType.Compliance, ProblemType.AppearanceWithMaxCompliance
                       ,ProblemType.ComplianceWithSymmetry
                       ,ProblemType.AppearanceWithMaxComplianceAndSymmetry)
        

