# PowerShell script to automate SSH login to mie server
param(
    [string]$Command = "cd WORK && pwd"
)

$password = "10428Kou"
$username = "nagasaku"
$hostname = "mie"

# Using plink (PuTTY) if available, otherwise standard ssh with echo
if (Get-Command plink -ErrorAction SilentlyContinue) {
    echo $password | plink -ssh -pw $password $username@$hostname $Command
} else {
    # For OpenSSH on Windows
    $securePassword = ConvertTo-SecureString $password -AsPlainText -Force
    $credential = New-Object System.Management.Automation.PSCredential ($username, $securePassword)

    # Execute command via SSH
    ssh $username@$hostname $Command
}