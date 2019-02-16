$system = Get-WmiObject Win32_ComputerSystem -EnableAllPrivileges
if ($system.AutomaticManagedPagefile) {
    Write-Host "Disabling managed page file settings"
    $system.AutomaticManagedPagefile = $false
    $system.Put() | Out-Null
}

$new_size = 3072
Write-Host "Changing page file size to $new_size MB"
$pagefile_settings = Get-WmiObject -Class Win32_PageFileSetting
$pagefile_settings.InitialSize = $new_size
$pagefile_settings.MaximumSize = $new_size
$pagefile_settings.Put() | Out-Null
